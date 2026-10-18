/* egaExplorer: search, gene view with arc view and neighbourhood graph, browse views. No network access. */
(function () {
  'use strict';

  var INDEX = JSON.parse(document.getElementById('ega-index').textContent);
  var D = null;          // the embedded data, after decompression
  var G = null;          // prepared graph
  var MINUS = '−';
  var SVGNS = 'http://www.w3.org/2000/svg';
  var state = { view: 'overview', gene: null, k: null, maxNodes: null, mode: 'change', mark: null, selected: null };
  var geneLower = {};
  INDEX.genes.forEach(function (g) { geneLower[g[0].toLowerCase()] = g[0]; });

  // ------------------------------------------------------------------------------------------------------------------
  // small helpers

  function $(sel, el) { return (el || document).querySelector(sel); }
  function esc(s) {
    return String(s == null ? '' : s).replace(/[&<>"']/g, function (c) {
      return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c];
    });
  }
  function el(tag, attrs, html) {
    var e = document.createElement(tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      if (attrs[k] === null || attrs[k] === undefined || attrs[k] === false) return;
      if (k === 'text') e.textContent = attrs[k]; else e.setAttribute(k, attrs[k]);
    });
    if (html !== undefined) e.innerHTML = html;
    return e;
  }
  function svg(tag, attrs) {
    var e = document.createElementNS(SVGNS, tag);
    if (attrs) Object.keys(attrs).forEach(function (k) {
      if (attrs[k] !== null && attrs[k] !== undefined) e.setAttribute(k, attrs[k]);
    });
    return e;
  }
  function isnum(v) { return typeof v === 'number' && isFinite(v); }
  function roundSig(v, d) {
    if (!isnum(v) || v === 0) return v;
    var p = d - 1 - Math.floor(Math.log10(Math.abs(v)));
    var f = Math.pow(10, p);
    return Math.round(v * f) / f;
  }
  function sig(v, d) {
    if (!isnum(v)) return '';
    d = d || 2;
    var r = roundSig(v, d);
    if (r === 0) return '0';
    var dec = Math.max(0, d - 1 - Math.floor(Math.log10(Math.abs(r))));
    return r.toLocaleString('en-US', { minimumFractionDigits: dec, maximumFractionDigits: dec }).replace('-', MINUS);
  }
  function log2txt(v) {
    if (!isnum(v)) return '';
    var r = Math.round(v * 10) / 10;
    if (r === 0) return '0.0';
    return (r > 0 ? '+' : MINUS) + Math.abs(r).toFixed(1);
  }
  function pq(v) {
    if (!isnum(v)) return '';
    if (v === 0) return '0';
    if (v < 0.001) {
      var parts = v.toExponential(1).split('e');
      var ex = parseInt(parts[1], 10);
      return parts[0] + ' × 10<sup>' + (ex < 0 ? MINUS : '') + Math.abs(ex) + '</sup>';
    }
    var t = sig(v, 2);
    return t.indexOf('.') >= 0 ? t.replace(/0+$/, '').replace(/\.$/, '') : t;
  }
  function int(v) { return isnum(v) ? Math.round(v).toLocaleString('en-US') : ''; }
  function distance(b) {
    b = Math.abs(b);
    if (b < 1000) return b.toFixed(0) + ' bp';
    if (b < 10000) return (b / 1000).toFixed(1) + ' kb';
    if (b < 1e6) return (b / 1000).toFixed(0) + ' kb';
    return (b / 1e6).toFixed(1) + ' Mb';
  }
  function joinWords(items, conj) {
    items = items.filter(Boolean);
    if (items.length <= 1) return items.join('');
    return items.slice(0, -1).join(', ') + ' ' + (conj || 'and') + ' ' + items[items.length - 1];
  }
  function capital(text) { return text.charAt(0).toUpperCase() + text.slice(1); }
  function geneTag(name) { return '<x-g>' + esc(name) + '</x-g>'; }
  function cssVar(name) { return getComputedStyle(document.documentElement).getPropertyValue(name).trim(); }
  function hexRgb(h) {
    h = h.replace('#', '');
    if (h.length === 3) h = h.split('').map(function (c) { return c + c; }).join('');
    return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)];
  }
  // diverging scale: lost blue, neutral grey midpoint, gained red (RdBu, colour-blind safe as a pair)
  var DIVERGING = null;
  function diverging(v) {
    if (!DIVERGING) {
      DIVERGING = [[-1, hexRgb('#2166ac')], [-0.5, hexRgb('#67a9cf')], [0, hexRgb(cssVar('--neutral') || '#f2f1ed')],
                   [0.5, hexRgb('#ef8a62')], [1, hexRgb('#b2182b')]];
    }
    if (!isnum(v)) return cssVar('--empty') || '#ffffff';
    var max = D.meta.maxLog2 || 2;
    var x = Math.max(-1, Math.min(1, v / max));
    for (var i = 1; i < DIVERGING.length; i++) {
      if (x <= DIVERGING[i][0]) {
        var a = DIVERGING[i - 1], b = DIVERGING[i];
        var t = (x - a[0]) / (b[0] - a[0]);
        var c = a[1].map(function (av, j) { return Math.round(av + t * (b[1][j] - av)); });
        return 'rgb(' + c.join(',') + ')';
      }
    }
    return 'rgb(178,24,43)';
  }

  // ------------------------------------------------------------------------------------------------------------------
  // tooltip

  var tip = document.getElementById('tip');
  function showTip(evt, html) {
    tip.innerHTML = html;
    tip.hidden = false;
    var pad = 12, w = tip.offsetWidth, h = tip.offsetHeight;
    var x = evt.clientX + pad, y = evt.clientY + pad;
    if (x + w > window.innerWidth - 8) x = Math.max(8, evt.clientX - w - pad);
    if (y + h > window.innerHeight - 8) y = Math.max(8, evt.clientY - h - pad);
    tip.style.left = x + 'px';
    tip.style.top = y + 'px';
  }
  function hideTip() { tip.hidden = true; }

  // ------------------------------------------------------------------------------------------------------------------
  // loading the embedded data

  function setStatus(text) {
    var s = document.getElementById('status');
    s.textContent = text || '';
    s.hidden = !text;
  }

  function load() {
    var text = document.getElementById('ega-data').textContent.trim();
    if (typeof DecompressionStream === 'undefined') {
      setStatus('This browser cannot decompress the embedded data (DecompressionStream is missing). Use a current ' +
                'version of Firefox, Chromium, Chrome, Edge or Safari.');
      return Promise.resolve();
    }
    var bin = atob(text);
    var bytes = new Uint8Array(bin.length);
    for (var i = 0; i < bin.length; i++) bytes[i] = bin.charCodeAt(i);
    var stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip'));
    return new Response(stream).text().then(function (json) {
      D = JSON.parse(json);
      prepare();
      setStatus('');
      route();
    }).catch(function (err) {
      setStatus('The embedded data could not be read: ' + err);
    });
  }

  function prepare() {
    var nodes = D.graph.nodes, edges = D.graph.edges;
    var n = nodes.s.length, m = edges.u.length;
    var deg = new Int32Array(n + 1);
    for (var e = 0; e < m; e++) { deg[edges.u[e] + 1]++; deg[edges.v[e] + 1]++; }
    for (var i = 0; i < n; i++) deg[i + 1] += deg[i];
    var fill = deg.slice(0, n);
    var adj = new Int32Array(2 * m);
    for (e = 0; e < m; e++) { adj[fill[edges.u[e]]++] = e; adj[fill[edges.v[e]]++] = e; }
    G = { n: n, m: m, off: deg, adj: adj, explains: {}, baitsByChrom: {} };
    var ignored = D.meta.ignored || [];
    G.ignored = ignored;
    // promoters per chromosome, for the nearest gene of an element
    for (i = 0; i < n; i++) {
      if (!nodes.b[i]) continue;
      var names = geneNames(nodes.g[i]);
      if (!names.length) continue;
      var c = nodes.c[i];
      (G.baitsByChrom[c] = G.baitsByChrom[c] || []).push([mid(i), names[0]]);
    }
    Object.keys(G.baitsByChrom).forEach(function (c) { G.baitsByChrom[c].sort(function (a, b) { return a[0] - b[0]; }); });
    // which genes each node explains (nodes that carry the explanation, egaReport)
    Object.keys(D.genes).forEach(function (gene) {
      (D.genes[gene].x || []).forEach(function (x) {
        (G.explains[x[0]] = G.explains[x[0]] || []).push([gene, x[1], x[2], x[3]]);
      });
    });
  }

  // ------------------------------------------------------------------------------------------------------------------
  // nodes and places

  function geneNames(text) {
    var names = String(text || '').split(',').filter(function (g) { return g && g !== 'nan'; });
    var informative = names.filter(function (g) {
      return !G.ignored.some(function (p) { return g.indexOf(p) === 0; });
    });
    return informative.length ? informative : names;
  }
  function mid(i) { return D.graph.nodes.s[i] + D.graph.nodes.e[i] / 2; }
  function chrom(i) { return D.graph.chroms[D.graph.nodes.c[i]]; }
  function coord(i) {
    var s = D.graph.nodes.s[i];
    return chrom(i) + ':' + s.toLocaleString('en-US') + '-' + (s + D.graph.nodes.e[i]).toLocaleString('en-US');
  }
  function classes(i) { return [D.graph.classes[D.graph.nodes.ca[i]], D.graph.classes[D.graph.nodes.cb[i]]]; }
  var ADJ = { 'active promoter': 'active', 'bivalent promoter': 'bivalent', 'poised promoter': 'poised',
              'primed promoter': 'primed', 'repressed promoter': 'repressed',
              'heterochromatin-like promoter': 'heterochromatin-like', 'architectural promoter': 'architectural',
              'unmarked promoter': 'unmarked' };
  function classNoun(c) { return !c || /marks of no role$/.test(c) ? 'element' : c; }
  function nearestGene(i, exclude) {
    var list = G.baitsByChrom[D.graph.nodes.c[i]];
    if (!list) return null;
    var x = mid(i), lo = 0, hi = list.length;
    while (lo < hi) { var h = (lo + hi) >> 1; if (list[h][0] < x) lo = h + 1; else hi = h; }
    var best = null;
    for (var k = Math.max(0, lo - 12); k < Math.min(list.length, lo + 12); k++) {
      if (list[k][1] === exclude) continue;
      var d = Math.abs(list[k][0] - x);
      if (!best || d < best[1]) best = [list[k][1], d];
    }
    return best;
  }
  /* A node named as in egaReport (Places.label): 'the bivalent promoter of Hoxa2' or 'an active enhancer 120 kb
     upstream of Hoxa3 (near Skap2)'. */
  function place(i, ref) {
    var nodes = D.graph.nodes, names = D.graph.nodes.b[i] ? geneNames(nodes.g[i]) : [];
    var cls = classes(i), tss = D.graph.tss[nodes.t[i]];
    var otherEnd = false;
    if (!names.length && tss === 'non-baited promoter') { names = geneNames(nodes.g[i]); otherEnd = names.length > 0; }
    if (names.length) {
      var suffix = otherEnd ? ' (an other end)' : '';
      if (cls[0] === cls[1] && ADJ[cls[0]]) return 'the ' + ADJ[cls[0]] + ' promoter of ' + geneTag(names[0]) + suffix;
      if (cls[0] === cls[1] || !cls[0]) return 'the ' + geneTag(names[0]) + ' promoter' + suffix;
      return 'the ' + geneTag(names[0]) + ' promoter (' + (otherEnd ? 'an other end; ' : '') +
        perCondition(cls, function (c) { return ADJ[c] || 'marks of no role'; }) + ')';
    }
    var noun = 'element', extra = '';
    if (cls[0] || cls[1]) {
      if (cls[0] === cls[1]) noun = classNoun(cls[0]); else extra = perCondition(cls, classNoun);
    }
    var near = nearestGene(i, ref ? ref.gene : null), text;
    if (ref && ref.node !== undefined && chrom(ref.node) === chrom(i)) {
      var dist = mid(i) - mid(ref.node), side = 'from';
      if (ref.strand === '+' || ref.strand === '-') side = ((dist < 0) === (ref.strand === '+')) ? 'upstream of' : 'downstream of';
      text = noun + ' ' + distance(dist) + ' ' + side + ' ' + geneTag(ref.gene);
      var inside = [];
      if (near && near[1] < Math.abs(dist)) inside.push('near ' + geneTag(near[0]));
      if (extra) inside.push(extra);
      if (inside.length) text += ' (' + inside.join('; ') + ')';
    } else {
      text = noun + ' on ' + esc(chrom(i)) + (near ? ' near ' + geneTag(near[0]) : '');
      if (extra) text += ' (' + extra + ')';
    }
    return (/^[aeiou]/i.test(text) ? 'an ' : 'a ') + text;
  }
  function perCondition(cls, words) {
    return cls.map(function (c, k) { return words(c) + ' in ' + esc(D.meta.names[k]); }).join(', ');
  }
  function markChange(i, mark) {
    var v = D.graph.nodes.m[mark];
    if (!v) return NaN;
    var a = v[0][i], b = v[1][i];
    if (!isnum(a) || !isnum(b)) return NaN;
    return Math.log2((b + 1) / (a + 1));
  }
  function edgeStrength(e, which) {
    var E = D.graph.edges, a = E.a[e] || 0, b = E.b[e] || 0;
    if (which === 'A') return a;
    if (which === 'B') return b;
    var s = E.s[e];
    return s === 1 ? b : s === 2 ? a : Math.max(a, b);
  }
  function edgeWidth(strength) {
    if (!strength || strength <= 0) return 0.5;
    return 0.5 + 3.4 * Math.max(0, Math.min(1, Math.log(strength) / Math.log(3)));
  }
  var STATUS = ['shared', 'gained', 'lost'];

  // ------------------------------------------------------------------------------------------------------------------
  // the neighbourhood of a gene (as egaReport.neighbourhood_subgraph)

  function neighbourhood(seeds, k, maxNodes, explaining) {
    var E = D.graph.edges, dist = new Map(), parent = new Map(), frontier = [];
    seeds.forEach(function (s) { dist.set(s, 0); frontier.push(s); });
    for (var step = 1; step <= Math.max(k, 1); step++) {
      var next = [];
      frontier.forEach(function (u) {
        for (var p = G.off[u]; p < G.off[u + 1]; p++) {
          var e = G.adj[p], w = E.u[e] === u ? E.v[e] : E.u[e];
          if (!dist.has(w)) { dist.set(w, step); parent.set(w, u); next.push(w); }
        }
      });
      frontier = next;
    }
    var keep, reduced = false;
    if (dist.size <= maxNodes) {
      keep = new Set(dist.keys());
    } else {
      reduced = true;
      keep = new Set(seeds);
      var strength = new Map();
      seeds.forEach(function (s) {
        for (var p = G.off[s]; p < G.off[s + 1]; p++) {
          var e = G.adj[p], w = E.u[e] === s ? E.v[e] : E.u[e];
          if (dist.get(w) === 1) strength.set(w, Math.max(strength.get(w) || 0, edgeStrength(e)));
        }
      });
      var partners = Array.from(strength.keys()).sort(function (a, b) { return strength.get(b) - strength.get(a) || a - b; });
      partners.slice(0, Math.max(10, Math.floor(maxNodes / 3))).forEach(function (w) { keep.add(w); });
      explaining.forEach(function (t) {
        var x = t;
        while (x !== undefined && !keep.has(x)) { keep.add(x); x = parent.get(x); }
      });
    }
    var c0 = new Set(seeds.map(function (s) { return D.graph.nodes.c[s]; }));
    keep.forEach(function (x) { if (!c0.has(D.graph.nodes.c[x])) keep.delete(x); });
    var edges = [], degree = new Map();
    keep.forEach(function (u) {
      for (var p = G.off[u]; p < G.off[u + 1]; p++) {
        var e = G.adj[p], w = E.u[e] === u ? E.v[e] : E.u[e];
        if (keep.has(w) && u < w) { edges.push(e); degree.set(u, (degree.get(u) || 0) + 1); degree.set(w, (degree.get(w) || 0) + 1); }
      }
    });
    var nodes = Array.from(keep).filter(function (x) { return degree.get(x) || seeds.indexOf(x) >= 0; });
    nodes.sort(function (a, b) { return mid(a) - mid(b); });
    return { nodes: nodes, edges: edges, dist: dist, reduced: reduced, total: dist.size };
  }

  // force-directed layout (Fruchterman and Reingold), deterministic start from the genomic order
  function layout(sub, width, height) {
    var nodes = sub.nodes, n = nodes.length, pos = {}, E = D.graph.edges;
    var seed = 7;
    function rand() { seed = (seed * 16807) % 2147483647; return seed / 2147483647; }
    nodes.forEach(function (x, i) {
      var a = 2 * Math.PI * i / Math.max(n, 1);
      pos[x] = [0.5 + 0.35 * Math.cos(a) + 0.02 * (rand() - 0.5), 0.5 + 0.35 * Math.sin(a) + 0.02 * (rand() - 0.5)];
    });
    if (n < 2) { if (n === 1) pos[nodes[0]] = [0.5, 0.5]; return pos; }
    var kk = Math.sqrt(1 / n) * 0.9, t = 0.1, iters = n > 200 ? 220 : 320;
    var disp = {};
    for (var it = 0; it < iters; it++) {
      nodes.forEach(function (v) { disp[v] = [0, 0]; });
      for (var i = 0; i < n; i++) {
        var v = nodes[i], pv = pos[v];
        for (var j = i + 1; j < n; j++) {
          var u = nodes[j], pu = pos[u];
          var dx = pv[0] - pu[0], dy = pv[1] - pu[1], d2 = dx * dx + dy * dy + 1e-9;
          var f = kk * kk / d2;
          disp[v][0] += dx * f; disp[v][1] += dy * f; disp[u][0] -= dx * f; disp[u][1] -= dy * f;
        }
      }
      sub.edges.forEach(function (e) {
        var a = E.u[e], b = E.v[e], pa = pos[a], pb = pos[b];
        var dx = pa[0] - pb[0], dy = pa[1] - pb[1], d = Math.sqrt(dx * dx + dy * dy) + 1e-9;
        var f = d / kk;
        disp[a][0] -= dx * f; disp[a][1] -= dy * f; disp[b][0] += dx * f; disp[b][1] += dy * f;
      });
      nodes.forEach(function (v) {
        var d = disp[v], len = Math.sqrt(d[0] * d[0] + d[1] * d[1]) + 1e-9, step = Math.min(len, t);
        pos[v][0] += d[0] / len * step + (0.5 - pos[v][0]) * 0.01;
        pos[v][1] += d[1] / len * step + (0.5 - pos[v][1]) * 0.01;
      });
      t = Math.max(0.002, t * 0.985);
    }
    // fit into the drawing area with the aspect of the panel
    var xs = nodes.map(function (v) { return pos[v][0]; }), ys = nodes.map(function (v) { return pos[v][1]; });
    var x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs), y0 = Math.min.apply(null, ys), y1 = Math.max.apply(null, ys);
    var out = {}, padX = 36, padY = 26;
    nodes.forEach(function (v) {
      out[v] = [padX + (pos[v][0] - x0) / ((x1 - x0) || 1) * (width - 2 * padX),
                padY + (pos[v][1] - y0) / ((y1 - y0) || 1) * (height - 2 * padY)];
    });
    // keep markers apart
    for (var r = 0; r < 60; r++) {
      var moved = false;
      for (i = 0; i < n; i++) for (j = i + 1; j < n; j++) {
        var p = out[nodes[i]], q = out[nodes[j]], ddx = q[0] - p[0], ddy = q[1] - p[1], dd = Math.sqrt(ddx * ddx + ddy * ddy);
        if (dd < 11) {
          var push = (11 - dd) / 2 + 0.1, ux = dd ? ddx / dd : 1, uy = dd ? ddy / dd : 0;
          p[0] -= ux * push; p[1] -= uy * push; q[0] += ux * push; q[1] += uy * push; moved = true;
        }
      }
      if (!moved) break;
    }
    nodes.forEach(function (v) {
      out[v][0] = Math.max(8, Math.min(width - 8, out[v][0]));
      out[v][1] = Math.max(8, Math.min(height - 8, out[v][1]));
    });
    return out;
  }

  // ------------------------------------------------------------------------------------------------------------------
  // routing

  function parseHash() {
    var out = {};
    location.hash.replace(/^#/, '').split('&').forEach(function (part) {
      if (!part) return;
      var kv = part.split('=');
      out[decodeURIComponent(kv[0])] = decodeURIComponent(kv.slice(1).join('=') || '');
    });
    return out;
  }
  function setHash(params, replace) {
    var text = '#' + Object.keys(params).filter(function (k) { return params[k] !== null && params[k] !== undefined; })
      .map(function (k) { return encodeURIComponent(k) + '=' + encodeURIComponent(params[k]); }).join('&');
    if (location.hash === text) { route(); return; }
    if (replace) { history.replaceState(null, '', text); route(); } else location.hash = text;
  }
  function route() {
    var p = parseHash();
    if (!D) return;
    hideTip();
    if (p.gene) {
      var name = D.genes[p.gene] ? p.gene : geneLower[p.gene.toLowerCase()];
      if (name) {
        var changed = state.gene !== name;
        state.gene = name;
        state.view = 'gene';
        if (changed) { state.selected = null; state.mark = null; }
        state.k = p.k ? Math.max(1, Math.min(3, parseInt(p.k, 10) || 1)) : null;
        if (p.node !== undefined && p.node !== '') state.selected = parseInt(p.node, 10);
        if (changed) document.getElementById('q').value = name;
        renderGene(changed);
        markTab();
        return;
      }
      state.view = 'overview';
      renderOverview('No gene named "' + p.gene + '" in this comparison.');
      markTab();
      return;
    }
    state.view = p.view || 'overview';
    if (state.view === 'gene' && state.gene) { renderGene(false); markTab(); return; }
    var views = { overview: renderOverview, genes: renderGenes, communities: renderCommunities,
                  switching: renderSwitching, triggers: renderTriggers };
    (views[state.view] || renderOverview)();
    markTab();
  }
  function markTab() {
    var tabGene = document.getElementById('tab-gene');
    tabGene.hidden = !state.gene;
    if (state.gene) tabGene.setAttribute('href', '#gene=' + encodeURIComponent(state.gene));
    document.querySelectorAll('#tabs a').forEach(function (a) {
      a.classList.toggle('on', a.getAttribute('data-view') === state.view);
      if (a.getAttribute('data-view') === state.view) a.setAttribute('aria-current', 'page'); else a.removeAttribute('aria-current');
    });
  }
  window.addEventListener('hashchange', route);

  // gene names in text are links to the gene
  function linkGenes(root) {
    root.querySelectorAll('x-g').forEach(function (g) {
      var name = g.textContent;
      if (D.genes[name] && name !== state.gene) {
        g.classList.add('link');
        g.setAttribute('tabindex', '0');
        g.setAttribute('role', 'link');
      }
    });
  }
  document.addEventListener('click', function (evt) {
    var g = evt.target.closest && evt.target.closest('x-g.link');
    if (g) { evt.preventDefault(); setHash({ gene: g.textContent }); }
  });
  document.addEventListener('keydown', function (evt) {
    if (evt.key === 'Enter' && evt.target.matches && evt.target.matches('x-g.link')) setHash({ gene: evt.target.textContent });
  });

  // ------------------------------------------------------------------------------------------------------------------
  // search with autocomplete

  (function search() {
    var input = document.getElementById('q'), list = document.getElementById('suggest');
    var items = [], active = -1;
    function close() { list.hidden = true; input.setAttribute('aria-expanded', 'false'); active = -1; }
    function matches(text) {
      var t = text.trim().toLowerCase();
      if (!t) return [];
      var starts = [], contains = [];
      for (var i = 0; i < INDEX.genes.length; i++) {
        var g = INDEX.genes[i], low = g[0].toLowerCase();
        if (low.indexOf(t) === 0) starts.push(g); else if (low.indexOf(t) > 0) contains.push(g);
        if (starts.length >= 40) break;
      }
      // exact name, then differential genes, then the rest
      starts.sort(function (a, b) { return (b[0].toLowerCase() === t) - (a[0].toLowerCase() === t) || b[1] - a[1] || a[0].length - b[0].length; });
      contains.sort(function (a, b) { return b[1] - a[1] || a[0].length - b[0].length; });
      return starts.concat(contains).slice(0, 12);
    }
    function show() {
      items = matches(input.value);
      list.innerHTML = '';
      if (!input.value.trim()) { close(); return; }
      if (!items.length) {
        list.appendChild(el('li', { class: 'none', role: 'option', text: 'No gene of this comparison matches.' }));
      }
      items.forEach(function (g, i) {
        var li = el('li', { role: 'option', id: 'sg' + i, 'aria-selected': i === active ? 'true' : 'false' });
        li.innerHTML = '<span><x-g>' + esc(g[0]) + '</x-g></span><span class="hint">' +
          (g[1] ? 'differential, log2 ' + log2txt(g[2]) : 'not differential') + '</span>';
        li.addEventListener('mousedown', function (evt) { evt.preventDefault(); pick(g[0]); });
        list.appendChild(li);
      });
      list.hidden = false;
      input.setAttribute('aria-expanded', 'true');
    }
    function pick(name) { close(); input.value = name; input.blur(); setHash({ gene: name }); }
    input.addEventListener('input', function () { active = -1; show(); });
    input.addEventListener('focus', function () { if (input.value.trim() && input.value !== state.gene) show(); });
    input.addEventListener('blur', function () { setTimeout(close, 120); });
    input.addEventListener('keydown', function (evt) {
      if (evt.key === 'ArrowDown' || evt.key === 'ArrowUp') {
        if (list.hidden) show();
        active = Math.max(0, Math.min(items.length - 1, active + (evt.key === 'ArrowDown' ? 1 : -1)));
        Array.prototype.forEach.call(list.children, function (li, i) { li.setAttribute('aria-selected', i === active ? 'true' : 'false'); });
        if (items.length) input.setAttribute('aria-activedescendant', 'sg' + active);
        evt.preventDefault();
      } else if (evt.key === 'Enter') {
        var name = active >= 0 && items[active] ? items[active][0] : (geneLower[input.value.trim().toLowerCase()] || (items[0] && items[0][0]));
        if (name) pick(name);
        evt.preventDefault();
      } else if (evt.key === 'Escape') { close(); }
    });
  })();

  // ------------------------------------------------------------------------------------------------------------------
  // generic sortable, filterable table

  function makeTable(host, opts) {
    var sortKey = opts.sort || null, sortDir = opts.sortDir || 1, limit = opts.page || 100, filtered = opts.rows;
    var wrap = el('div', { class: 'table-scroll' }), table = el('table'), foot = el('div', { class: 'table-foot' });
    var count = el('span', { class: 'count' }), more = el('button', { class: 'plain', type: 'button', text: 'Show more rows' });
    wrap.appendChild(table);
    host.appendChild(wrap);
    foot.appendChild(count);
    foot.appendChild(more);
    host.appendChild(foot);
    more.addEventListener('click', function () { limit += opts.page || 100; draw(); });
    function value(col, row) { return col.value ? col.value(row) : row[col.key]; }
    function draw() {
      var rows = filtered.slice();
      if (sortKey !== null) {
        var col = opts.columns.filter(function (c) { return c.id === sortKey; })[0];
        rows.sort(function (a, b) {
          var x = value(col, a), y = value(col, b);
          var nx = x === null || x === undefined || x === '' || (typeof x === 'number' && !isFinite(x));
          var ny = y === null || y === undefined || y === '' || (typeof y === 'number' && !isFinite(y));
          if (nx || ny) return nx - ny;
          return (x < y ? -1 : x > y ? 1 : 0) * sortDir;
        });
      }
      var head = '<thead><tr>' + opts.columns.map(function (c) {
        var arrow = c.id === sortKey ? '<span class="arrow">' + (sortDir > 0 ? '▲' : '▼') + '</span>' : '';
        return '<th class="sort' + (c.num ? ' num' : '') + '" data-id="' + c.id + '"' + (c.title ? ' title="' + esc(c.title) + '"' : '') +
          ' scope="col" aria-sort="' + (c.id === sortKey ? (sortDir > 0 ? 'ascending' : 'descending') : 'none') + '">' + esc(c.label) + arrow + '</th>';
      }).join('') + '</tr></thead>';
      var body = rows.slice(0, limit).map(function (r, k) {
        return '<tr data-k="' + k + '"' + (opts.onRow ? ' class="click"' : '') + '>' + opts.columns.map(function (c) {
          var v = c.render ? c.render(r) : esc(value(c, r));
          return '<td class="' + (c.num ? 'num' : '') + (c.cls ? ' ' + c.cls : '') + '">' + (v === undefined ? '' : v) + '</td>';
        }).join('') + '</tr>';
      }).join('');
      table.innerHTML = head + '<tbody>' + body + '</tbody>';
      table._rows = rows;
      count.textContent = (rows.length > limit ? 'Showing ' + int(limit) + ' of ' : '') + int(rows.length) + ' ' +
        (rows.length === 1 ? opts.noun[0] : opts.noun[1]) + (rows.length !== opts.rows.length ? ' (of ' + int(opts.rows.length) + ')' : '');
      more.hidden = rows.length <= limit;
      linkGenes(table);
    }
    table.addEventListener('click', function (evt) {
      var th = evt.target.closest('th.sort');
      if (th) {
        var id = th.getAttribute('data-id');
        if (sortKey === id) sortDir = -sortDir; else { sortKey = id; sortDir = 1; }
        draw();
        return;
      }
      if (evt.target.closest('x-g.link')) return;
      var tr = evt.target.closest('tr[data-k]');
      if (tr && opts.onRow) opts.onRow(table._rows[parseInt(tr.getAttribute('data-k'), 10)], tr);
    });
    draw();
    return { filter: function (fn) { filtered = opts.rows.filter(fn); limit = opts.page || 100; draw(); } };
  }

  function controlsRow(host) { var c = el('div', { class: 'controls' }); host.appendChild(c); return c; }
  function selectControl(host, label, options, value, onChange) {
    var lab = el('label'), sel = el('select');
    lab.appendChild(document.createTextNode(label));
    options.forEach(function (o) {
      var opt = el('option', { value: o[0], text: o[1] });
      if (String(o[0]) === String(value)) opt.selected = true;
      sel.appendChild(opt);
    });
    sel.addEventListener('change', function () { onChange(sel.value); });
    lab.appendChild(sel);
    host.appendChild(lab);
    return sel;
  }
  function checkControl(host, label, onChange) {
    var lab = el('label'), box = el('input', { type: 'checkbox' });
    lab.appendChild(box);
    lab.appendChild(document.createTextNode(label));
    box.addEventListener('change', function () { onChange(box.checked); });
    host.appendChild(lab);
    return box;
  }
  function textControl(host, label, placeholder, onChange) {
    var lab = el('label'), box = el('input', { type: 'search', placeholder: placeholder, autocomplete: 'off' });
    lab.appendChild(document.createTextNode(label));
    lab.appendChild(box);
    box.addEventListener('input', function () { onChange(box.value.trim().toLowerCase()); });
    host.appendChild(lab);
    return box;
  }
  function segControl(host, label, options, value, onChange) {
    var lab = el('label'), seg = el('span', { class: 'seg', role: 'group', 'aria-label': label });
    lab.appendChild(document.createTextNode(label));
    options.forEach(function (o) {
      var b = el('button', { type: 'button', 'aria-pressed': String(o[0]) === String(value) ? 'true' : 'false', text: o[1] });
      b.addEventListener('click', function () { onChange(o[0]); });
      seg.appendChild(b);
    });
    lab.appendChild(seg);
    host.appendChild(lab);
  }
  function view() { var v = document.getElementById('view'); v.innerHTML = ''; window.scrollTo(0, 0); return v; }
  function section(host, title, id) { var h = el('h3', { id: id, text: title }); host.appendChild(h); return h; }
  function notesFor(block) {
    return (block && block.n || []).map(function (k) { return '<p class="note">' + D.meta.notes[k] + '</p>'; }).join('');
  }
  function levelText(lv) { return D.meta.levelLabels[lv] || lv; }
  function levelShort(lv) { return D.meta.levelShort[lv] || lv; }

  // ------------------------------------------------------------------------------------------------------------------
  // overview

  function exampleGene() {
    var best = null;
    INDEX.genes.forEach(function (g) { if (g[1] && isnum(g[2]) && (!best || Math.abs(g[2]) > Math.abs(best[2]))) best = g; });
    return best ? best[0] : (INDEX.genes.length ? INDEX.genes[0][0] : 'NAME');
  }

  function renderOverview(message) {
    var v = view(), M = D.meta;
    if (typeof message === 'string') v.appendChild(el('p', { class: 'note', text: message }));
    v.appendChild(el('h2', { text: 'Overview' }));
    var stats = el('div', { class: 'stats' });
    [[int(M.nDifferential), 'differential promoter fragments'], [int(M.nGenes), 'genes in the search index'],
     [int(M.nNodes), 'fragments in the graph'], [int(M.nEdges), 'interactions'],
     [int(M.edgeCounts.gained || 0), 'interactions gained in ' + M.names[1]],
     [int(M.edgeCounts.lost || 0), 'interactions lost in ' + M.names[1]]].forEach(function (s) {
      stats.appendChild(el('div', null, '<div class="v">' + s[0] + '</div><div class="k">' + esc(s[1]) + '</div>'));
    });
    v.appendChild(stats);
    var cascade = M.cascade || {};
    if (Object.keys(cascade).length) {
      section(v, 'Nearest level that explains a gene');
      var host = el('div');
      v.appendChild(host);
      var main = function (k) { return k === 'gene' || k === 'unexplained' || /^step\d+$/.test(k); };
      var total = Object.keys(cascade).filter(main).reduce(function (s, k) { return s + cascade[k]; }, 0);
      var label = function (k) { return k === 'gene' ? 'Promoter' : k === 'step1' ? 'Direct partners' : k === 'unexplained' ? 'Not explained' : k.replace('step', '') + ' steps away'; };
      var parts = { marks_only: 'histone marks only', interactions_only: 'own contacts only', both: 'histone marks and own contacts',
                    state: 'chromatin state', poised_unresolved: 'poised, trigger not measured', marks: 'histone marks' };
      var rows = Object.keys(cascade).map(function (k) {
        if (main(k)) return [label(k), cascade[k], cascade[k] / total, false];
        var m = /^(gene|step\d+)_(.+)$/.exec(k);
        return [m ? 'of which: ' + (parts[m[2]] || m[2].replace(/_/g, ' ')) : k, cascade[k], cascade[k] / total, true];
      });
      var html = '<div class="table-scroll"><table><thead><tr><th scope="col">Level</th><th class="num" scope="col">Genes</th><th class="num" scope="col">Fraction</th></tr></thead><tbody>' +
        rows.map(function (r) { return '<tr' + (r[3] ? ' class="muted"' : '') + '><td' + (r[3] ? ' style="padding-left:1.4rem"' : '') + '>' + esc(r[0]) + '</td><td class="num">' + int(r[1]) + '</td><td class="num">' + (100 * r[2]).toFixed(1) + '%</td></tr>'; }).join('') +
        '</tbody></table></div><p class="caption"><span class="label">Table 1.</span> Differential genes (' + int(total) + ') by the nearest level that explains their expression change (explanation cascade of egaDifferential), with the percentage of all differential genes; indented rows split a level by the kind of explanation.</p>';
      host.innerHTML = html;
    }
    section(v, 'Using the explorer');
    v.appendChild(el('ul', { class: 'plain' },
      '<li>Search a gene in the box at the top left; suggestions list differential genes first. The address of a gene, ' +
      'for example <code>#gene=' + esc(exampleGene()) + '</code>, opens it directly and can be shared with the file.</li>' +
      '<li>A gene page holds the text of the gene report (egaReport): at a glance, the explanation ladder, the evidence by level, chromatin states, element classes, hub and community.</li>' +
      '<li>The arc view and the neighbourhood graph are drawn from the embedded graph for 1 to 3 steps. Hover a node or an interaction for its details; click a node to see what it explains.</li>' +
      '<li>The browse views list the differential genes, the communities and their transitions' +
      (D.switching ? ', the enhancer switching events' : '') + (D.triggers ? ' and the candidate triggers' : '') +
      '; their columns sort on click and gene names open the gene.</li>'));
    section(v, 'Comparison');
    v.appendChild(el('p', null, 'Condition A: ' + esc(M.labels[0]) + (M.labels[0] !== M.names[0] ? ' (' + esc(M.names[0]) + ')' : '') +
      '; condition B: ' + esc(M.labels[1]) + (M.labels[1] !== M.names[1] ? ' (' + esc(M.names[1]) + ')' : '') +
      '. Changes are log2 fold changes ' + esc(M.names[1]) + ' / ' + esc(M.names[0]) + '. Marks and declared roles: ' +
      M.marks.map(function (m) { return esc(m) + ' ' + esc(M.roles[m] || ''); }).join(', ') +
      '. Differential genes: |log2 fold change| of at least ' + sig(M.rnaThreshold, 2) + '; a level explains a gene when its change agrees with the expression change by at least log2 ' +
      sig(M.chromatinThreshold, 2) + '. Levels from ' + esc(M.levelsSource) + '.'));
  }

  // ------------------------------------------------------------------------------------------------------------------
  // gene page

  function geneSeeds(gene) { return D.genes[gene].n || []; }

  function renderGene(changed) {
    var gene = state.gene, info = D.genes[gene], M = D.meta;
    var v = changed || !document.getElementById('gene-page') ? view() : document.getElementById('view');
    if (!changed && document.getElementById('gene-page')) { drawFigures(); return; }
    var page = el('div', { id: 'gene-page' });
    v.appendChild(page);
    var title = el('h2', null, '<x-g>' + esc(gene) + '</x-g><span class="pill">' + (info.d ? 'differential' : 'not differential') + '</span>');
    page.appendChild(title);
    var seeds = geneSeeds(gene);
    page.appendChild(el('p', { class: 'sub' }, (seeds.length === 1 ? 'Promoter fragment ' : 'Promoter fragments ') +
      seeds.map(function (s) { return '<span class="mono">' + esc(coord(s)) + '</span>'; }).join(', ') +
      (isnum(info.r) ? '; expression log2 ' + log2txt(info.r) + ' (' + esc(M.names[1]) + ' / ' + esc(M.names[0]) + ')' : '')));

    section(page, 'At a glance', 'glance');
    var grid = el('div', { class: 'glance-grid' + (info.ladder ? '' : ' single') });
    grid.appendChild(el('p', { class: 'glance' }, info.glance));
    if (info.ladder) grid.appendChild(ladderFigure(info.ladder, gene));
    page.appendChild(grid);

    section(page, 'Interactions', 'interactions');
    var controls = controlsRow(page);
    controls.id = 'graph-controls';
    var arcs = el('figure', { id: 'arc-figure' });
    page.appendChild(arcs);
    var gfig = el('figure', { id: 'graph-figure' });
    page.appendChild(gfig);

    if (info.evidence && info.evidence.length) {
      section(page, 'Evidence by level', 'evidence');
      info.evidence.forEach(function (frag, k) {
        if (frag.h) page.appendChild(el('h4', null, frag.h));
        var html = '<div class="table-scroll"><table><thead><tr>' + frag.cols.map(function (c) {
          return '<th scope="col"' + (c === 'P' || c === 'q' ? ' class="num"' : '') + '>' + esc(c) + '</th>';
        }).join('') + '</tr></thead><tbody>' + frag.rows.map(function (r, i) {
          return '<tr class="' + esc(frag.cls[i] || '') + '">' + r.map(function (c, j) {
            var col = frag.cols[j];
            var cls = col === 'What changed' ? 'long' : (col === 'P' || col === 'q') ? 'num' : (col === 'Level' || col === 'Explains' || col === 'Distance control') ? 'key' : '';
            return '<td class="' + cls + '">' + c + '</td>';
          }).join('') + '</tr>';
        }).join('') + '</tbody></table></div>';
        page.appendChild(el('div', null, html));
        if (k === 0) {
          page.appendChild(el('p', { class: 'caption' }, '<span class="label">Evidence.</span> From the promoter outwards. Mark changes are log2 fold changes ' +
            esc(M.names[1]) + ' / ' + esc(M.names[0]) + ', at the promoter itself or over the nodes of a step; explains: the change agrees with the expression change by at least log2 ' +
            sig(M.chromatinThreshold, 2) + '. P and q: per-gene permutation test of the step; distance control: the change exceeds that at distance-matched mirror loci. ' +
            'A blue bar marks a level that explains the gene, is significant and passes the distance control.'));
        }
        if (frag.learned && frag.learned.length) {
          page.appendChild(el('p', { class: 'note' }, 'Under the mark roles learned from these data, ' + frag.learned.join('; ') + '.'));
        }
      });
    }
    if (info.tads) {
      section(page, 'Domain boundaries', 'boundaries');
      page.appendChild(el('div', null, info.tads.h + notesFor(info.tads)));
    }
    if (info.states) {
      section(page, 'Chromatin states', 'states');
      page.appendChild(el('div', null, info.states.h + notesFor(info.states)));
    }
    if (info.elements) {
      section(page, 'Element classes', 'elements');
      page.appendChild(el('div', null, info.elements.h.replace(/^<h4>[^<]*<\/h4>/, '') + notesFor(info.elements)));
    }
    if (info.hub || info.community) {
      section(page, 'Hub and community', 'hub');
      if (info.hub) page.appendChild(el('div', null, info.hub.h));
      if (info.community) page.appendChild(el('div', null, info.community.h));
    }
    geneTriggers(page, gene);
    geneSwitching(page, gene);
    linkGenes(page);
    drawFigures();
  }

  function ladderFigure(lad, gene) {
    var fig = el('figure', { class: 'ladder' });
    var nCols = lad.cols.length, hasSupport = Object.keys(lad.support || {}).length > 0;
    var labelW = 236, colW = nCols > 1 ? 112 : 46, supportW = hasSupport ? 74 : 0, rowH = 24;
    var top = nCols > 1 || hasSupport ? 30 : 10;
    var W = labelW + nCols * colW + supportW + 8, H = top + lad.levels.length * rowH + 8;
    var s = svg('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, role: 'img', 'aria-label': 'Explanation ladder of ' + gene });
    s.style.width = '100%';
    s.style.maxWidth = W + 'px';
    lad.levels.forEach(function (lv, r) {
      var y = top + r * rowH + rowH / 2;
      var t = svg('text', { x: labelW - 10, y: y + 4, 'text-anchor': 'end', 'font-size': 12.5 });
      t.textContent = lad.labels[r];
      s.appendChild(t);
      if (hasSupport && lad.support[lv]) {
        var st = svg('text', { x: labelW + nCols * colW + 10, y: y + 4, 'font-size': 12 });
        st.textContent = lad.support[lv];
        s.appendChild(st);
      }
    });
    if (hasSupport) {
      var hd = svg('text', { x: labelW + nCols * colW + 10, y: 18, 'font-size': 11, class: 'muted' });
      hd.textContent = '3D support';
      s.appendChild(hd);
    }
    lad.cols.forEach(function (col, j) {
      var x = labelW + j * colW + colW / 2;
      if (nCols > 1) {
        var ct = svg('text', { x: x, y: 18, 'text-anchor': 'middle', 'font-size': 10.5, class: 'muted' });
        ct.textContent = col.l.length > 15 ? col.l.slice(0, 14) + '…' : col.l;
        var title = svg('title');
        title.textContent = col.l;
        ct.appendChild(title);
        s.appendChild(ct);
      }
      var ys = lad.levels.map(function (lv, r) { return col.c[lv] ? top + r * rowH + rowH / 2 : null; }).filter(function (y) { return y !== null; });
      if (ys.length > 1) s.appendChild(svg('line', { x1: x, x2: x, y1: Math.min.apply(null, ys), y2: Math.max.apply(null, ys), class: 'lad-line' }));
      lad.levels.forEach(function (lv, r) {
        var cell = col.c[lv];
        if (!cell) return;
        var y = top + r * rowH + rowH / 2;
        if (cell[1]) s.appendChild(svg('circle', { cx: x, cy: y, r: 9, class: 'lad-ring' }));
        s.appendChild(svg('circle', { cx: x, cy: y, r: 5.2, class: cell[0] ? 'lad-yes' : 'lad-no' }));
      });
    });
    fig.appendChild(s);
    fig.appendChild(el('figcaption', null, '<span class="label">Explanation ladder.</span> One row per level from the promoter outwards' +
      (nCols > 1 ? ', one column per promoter fragment' : '') + '. Filled: the level explains the expression change; open: it does not; ' +
      'blue ring: significant (q &lt; 0.05) and beyond distance-matched mirror loci.' +
      (hasSupport ? ' 3D support: how much closer graph neighbours at this step are in the nucleus than distance-matched pairs (0.5: no relation).' : '')));
    return fig;
  }

  function explainingFor(gene, k) {
    return (D.genes[gene].x || []).filter(function (x) { return x[1] <= k; });
  }
  function defaultK(gene) {
    var steps = (D.genes[gene].x || []).map(function (x) { return x[1]; });
    return Math.max(1, Math.min(3, steps.length ? Math.max.apply(null, steps) : 1));
  }

  function drawFigures() {
    var gene = state.gene, info = D.genes[gene], M = D.meta;
    var k = state.k || defaultK(gene), maxNodes = state.maxNodes || M.maxNodes;
    var mark = state.mark || info.mark || M.marks[0];
    var controls = document.getElementById('graph-controls');
    controls.innerHTML = '';
    segControl(controls, 'Steps', [[1, '1'], [2, '2'], [3, '3']], k, function (x) {
      setHash({ gene: gene, k: x }, true);
    });
    selectControl(controls, 'Largest drawn', [[40, '40 nodes'], [80, '80 nodes'], [150, '150 nodes'], [300, '300 nodes']],
      maxNodes, function (x) { state.maxNodes = parseInt(x, 10); drawFigures(); });
    selectControl(controls, 'Graph interactions', [['change', 'change'], ['A', M.names[0]], ['B', M.names[1]]], state.mode,
      function (x) { state.mode = x; drawFigures(); });
    if (M.marks.length) {
      selectControl(controls, 'Other ends coloured by', M.marks.map(function (m) { return [m, m]; }), mark,
        function (x) { state.mark = x; drawFigures(); });
    }
    var seeds = geneSeeds(gene).filter(function (s) { return s < G.n; });
    var arcsHost = document.getElementById('arc-figure'), graphHost = document.getElementById('graph-figure');
    arcsHost.innerHTML = '';
    graphHost.innerHTML = '';
    if (!seeds.length) { arcsHost.appendChild(el('p', { class: 'note', text: 'The promoter of this gene is not in the graph.' })); return; }
    var expl = explainingFor(gene, k);
    var explNodes = expl.map(function (x) { return x[0]; });
    var sub = neighbourhood(seeds, k, maxNodes, explNodes);
    var inSub = new Set(sub.nodes);
    expl = expl.filter(function (x) { return inSub.has(x[0]); });
    var ctx = { gene: gene, seeds: seeds, sub: sub, k: k, mark: mark, expl: expl,
                explSet: new Set(expl.map(function (x) { return x[0]; })), ref: { node: seeds[0], gene: gene, strand: info.strand } };
    if (!sub.edges.length) {
      arcsHost.appendChild(el('p', { class: 'note' }, '<x-g>' + esc(gene) + '</x-g> has no interactions in either ' + esc(M.kind === 'condition' ? 'condition' : M.kind) + '.'));
      return;
    }
    drawArcs(arcsHost, ctx);
    drawGraph(graphHost, ctx);
  }

  // node tooltip and node panel ----------------------------------------------------------------------------------
  function nodeSummary(i, ctx) {
    var N = D.graph.nodes, M = D.meta, cls = classes(i);
    var lines = ['<div class="t">' + capital(place(i, ctx.ref)) + '</div>', '<div class="m mono">' + esc(coord(i)) + '</div>'];
    var d = ctx.sub.dist.get(i);
    if (d !== undefined) lines.push('<div>' + (d === 0 ? 'promoter of ' + geneTag(ctx.gene) : d + (d === 1 ? ' step' : ' steps') + ' from the promoter') + '</div>');
    if (N.g[i]) lines.push('<div>Genes: ' + geneNames(N.g[i]).slice(0, 6).map(geneTag).join(', ') + '</div>');
    if (cls[0] || cls[1]) lines.push('<div>' + (cls[0] === cls[1] ? esc(cls[0]) + ' in both' : esc(cls[0] || 'n.a.') + ' in ' + esc(M.names[0]) + ', ' + esc(cls[1] || 'n.a.') + ' in ' + esc(M.names[1])) + '</div>');
    if (N.b[i] && isnum(N.r[i])) lines.push('<div>Expression log2 ' + log2txt(N.r[i]) + '</div>');
    var marks = M.marks.map(function (m) { var c = markChange(i, m); return isnum(c) ? esc(m) + ' ' + log2txt(c) : ''; }).filter(Boolean);
    if (marks.length) lines.push('<div>Marks, log2: ' + marks.join(', ') + '</div>');
    var x = ctx.expl.filter(function (e) { return e[0] === i; })[0];
    if (x) lines.push('<div>Carries the explanation ' + (x[1] === 1 ? 'at the direct partners' : x[1] + ' steps away') + ' (' + esc(x[2]) + ' ' + log2txt(x[3]) + ')</div>');
    return lines.join('');
  }

  function nodePanel(host, i, ctx) {
    var N = D.graph.nodes, M = D.meta;
    host.innerHTML = '';
    if (i === null || i === undefined || !ctx.sub.dist.has(i)) {
      host.appendChild(el('p', { class: 'note', text: 'Click a node to show its details and what it explains.' }));
      return;
    }
    var html = '<h4>' + capital(place(i, ctx.ref)) + '</h4><p class="mono">' + esc(coord(i)) + ' (' + distance(N.e[i]) + ')</p>';
    var d = ctx.sub.dist.get(i);
    html += '<p>' + (d === 0 ? 'Promoter of ' + geneTag(ctx.gene) + '.' : d + (d === 1 ? ' step' : ' steps') + ' from the promoter of ' + geneTag(ctx.gene) + '.') + '</p>';
    var cls = classes(i);
    if (cls[0] || cls[1]) html += '<p>Element class: ' + esc(cls[0] || 'n.a.') + ' in ' + esc(M.names[0]) + ', ' + esc(cls[1] || 'n.a.') + ' in ' + esc(M.names[1]) + '.</p>';
    if (N.b[i] && isnum(N.r[i])) html += '<p>Expression log2 ' + log2txt(N.r[i]) + '.</p>';
    if (M.marks.length) {
      html += '<div class="table-scroll"><table><thead><tr><th scope="col">Mark</th><th class="num" scope="col">' + esc(M.names[0]) +
        '</th><th class="num" scope="col">' + esc(M.names[1]) + '</th><th class="num" scope="col">log2</th></tr></thead><tbody>' +
        M.marks.map(function (m) {
          var v = N.m[m];
          return '<tr><td>' + esc(m) + '</td><td class="num">' + sig(v[0][i], 2) + '</td><td class="num">' + sig(v[1][i], 2) + '</td><td class="num">' + log2txt(markChange(i, m)) + '</td></tr>';
        }).join('') + '</tbody></table></div><p class="note">Fold enrichment over input; log2 of (' + esc(M.names[1]) + ' + 1) / (' + esc(M.names[0]) + ' + 1).</p>';
    }
    var counts = [0, 0, 0], E = D.graph.edges;
    for (var p = G.off[i]; p < G.off[i + 1]; p++) counts[E.s[G.adj[p]]]++;
    html += '<p>Interactions: ' + int(counts[0] + counts[2]) + ' in ' + esc(M.names[0]) + ', ' + int(counts[0] + counts[1]) + ' in ' + esc(M.names[1]) +
      ' (' + int(counts[1]) + ' gained, ' + int(counts[2]) + ' lost).</p>';
    var mine = ctx.expl.filter(function (e) { return e[0] === i; })[0];
    if (mine) html += '<p><strong>Explains ' + geneTag(ctx.gene) + '</strong> ' + (mine[1] === 1 ? 'at the direct partners' : mine[1] + ' steps away') + ': ' + esc(mine[2]) + ' ' + log2txt(mine[3]) + '.</p>';
    var others = (G.explains[i] || []).filter(function (e) { return e[0] !== ctx.gene; });
    if (others.length) {
      var shown = others.slice(0, 20).map(function (e) { return geneTag(e[0]) + ' (' + (e[1] === 1 ? 'step 1' : 'step ' + e[1]) + ', ' + esc(e[2]) + ' ' + log2txt(e[3]) + ')'; });
      html += '<p>Carries the explanation of ' + (others.length === 1 ? 'one other gene' : int(others.length) + ' other genes') + ': ' + shown.join(', ') + (others.length > 20 ? ' and ' + (others.length - 20) + ' more' : '') + '.</p>';
    } else if (!mine) {
      html += '<p class="note">This node carries no explanation of a differential gene.</p>';
    }
    var names = N.b[i] ? geneNames(N.g[i]).filter(function (g) { return D.genes[g] && g !== ctx.gene; }) : [];
    if (names.length) html += '<p>Open: ' + names.slice(0, 6).map(geneTag).join(', ') + '</p>';
    host.innerHTML = html;
    linkGenes(host);
  }

  // arc view ----------------------------------------------------------------------------------------------------------
  function drawArcs(host, ctx) {
    var M = D.meta, E = D.graph.edges, sub = ctx.sub;
    var W = Math.max(300, Math.min(host.clientWidth || 800, 1180));
    var narrow = W < 560;
    var hA = narrow ? 56 : 70, hC = narrow ? 64 : 84, gap = 22, labelH = 74, axisH = 40;
    var top = 16;
    var yA = top + hA, yB = yA + gap + hA, yC = yB + gap + 14 + hC;
    var H = yC + hC + labelH + axisH;
    var xs = sub.nodes.map(mid), x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
    var pad = 0.03 * ((x1 - x0) || 1), L = 12, R = W - 12;
    function X(v) { return L + (v - (x0 - pad)) / ((x1 + pad) - (x0 - pad) || 1) * (R - L); }
    var longest = 1;
    sub.edges.forEach(function (e) { longest = Math.max(longest, Math.abs(mid(E.u[e]) - mid(E.v[e]))); });
    var s = svg('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img',
                         'aria-label': 'Arc view of the interactions around ' + ctx.gene });
    function arc(e, y, h, up, cls, width) {
      var a = X(mid(E.u[e])), b = X(mid(E.v[e])), lo = Math.min(a, b), hi = Math.max(a, b);
      var rx = Math.max((hi - lo) / 2, 0.5), ry = Math.max(3, h * Math.sqrt(Math.abs(mid(E.u[e]) - mid(E.v[e])) / longest));
      var d = 'M' + lo.toFixed(1) + ' ' + y + ' A' + rx.toFixed(1) + ' ' + ry.toFixed(1) + ' 0 0 ' + (up ? 1 : 0) + ' ' + hi.toFixed(1) + ' ' + y;
      var hit = svg('path', { d: d, class: 'hit' });
      var p = svg('path', { d: d, class: 'arc ' + cls, 'stroke-width': width.toFixed(2) });
      hit.addEventListener('mousemove', function (evt) { showTip(evt, edgeTip(e, ctx)); });
      hit.addEventListener('mouseleave', hideTip);
      s.appendChild(hit);
      s.appendChild(p);
    }
    function axis(y) { s.appendChild(svg('line', { x1: L, x2: R, y1: y, y2: y, class: 'axis' })); }
    function label(x, y, text, cls) { var t = svg('text', { x: x, y: y, 'font-size': 11.5, class: cls || '' }); t.textContent = text; s.appendChild(t); }
    var countA = 0, countB = 0;
    sub.edges.forEach(function (e) { if (E.s[e] !== 1) countA++; if (E.s[e] !== 2) countB++; });
    label(L, top - 3, M.names[0] + ': ' + countA + ' interactions');
    sub.edges.forEach(function (e) { if (E.s[e] !== 1) arc(e, yA, hA - 6, true, 'e-cond', edgeWidth(edgeStrength(e, 'A'))); });
    axis(yA);
    label(L, yA + gap - 3, M.names[1] + ': ' + countB + ' interactions');
    sub.edges.forEach(function (e) { if (E.s[e] !== 2) arc(e, yB, hA - 6, true, 'e-cond', edgeWidth(edgeStrength(e, 'B'))); });
    axis(yB);
    var counts = [0, 0, 0];
    sub.edges.forEach(function (e) { counts[E.s[e]]++; });
    label(L, yB + gap + 8, 'Change: ' + counts[1] + ' gained in ' + M.names[1] + ' (above), ' + counts[2] + ' lost (below), ' + counts[0] + ' shared');
    [0, 1, 2].forEach(function (st) {
      sub.edges.forEach(function (e) { if (E.s[e] === st) arc(e, yC, hC - 4, st !== 2, 'e-' + STATUS[st], edgeWidth(edgeStrength(e))); });
    });
    axis(yC);
    // TAD boundaries of the span (egaExplorer --tads): dashed lines, thick for the strongest decile
    if (D.tads && sub.nodes.length) {
      var tadChrom = D.graph.chroms[D.graph.nodes.c[sub.nodes[0]]];
      (D.tads.chroms[tadChrom] || []).forEach(function (bd) {
        if (bd[0] < x0 - pad || bd[0] > x1 + pad) return;
        var xb = X(bd[0]);
        var line = svg('line', { x1: xb, x2: xb, y1: top, y2: yC + hC, class: 'tad-boundary' + (bd[2] ? ' strong' : '') });
        var title = svg('title', {});
        title.textContent = 'TAD boundary ' + tadChrom + ':' + (bd[0] / 1e6).toFixed(3) + ' Mb' + (bd[1] ? ', strength rank ' + bd[1] + ' of ' + D.tads.n : '');
        line.appendChild(title);
        s.appendChild(line);
      });
    }
    // markers on the axes: promoters, nodes that carry the explanation, the gene
    [yA, yB, yC].forEach(function (y) {
      sub.nodes.forEach(function (i) {
        var x = X(mid(i)), g;
        if (ctx.seeds.indexOf(i) >= 0) g = svg('rect', { x: x - 4.5, y: y - 4.5, width: 9, height: 9, class: 'seed-mark' });
        else if (ctx.explSet.has(i)) g = svg('circle', { cx: x, cy: y, r: 4, class: 'expl-mark' });
        else if (D.graph.nodes.b[i]) g = svg('rect', { x: x - 2.5, y: y - 2.5, width: 5, height: 5, class: 'bait-mark' });
        else return;
        g.style.cursor = 'pointer';
        g.addEventListener('mousemove', function (evt) { showTip(evt, nodeSummary(i, ctx)); });
        g.addEventListener('mouseleave', hideTip);
        g.addEventListener('click', function () { selectNode(i, ctx); });
        s.appendChild(g);
      });
    });
    // gene labels below the change, each gene once, the gene of the page first
    var placed = [], shown = new Set(), minGap = 12;
    var baits = sub.nodes.filter(function (i) { return D.graph.nodes.b[i]; });
    baits.sort(function (a, b) { return (ctx.seeds.indexOf(b) >= 0) - (ctx.seeds.indexOf(a) >= 0) || (ctx.explSet.has(b) - ctx.explSet.has(a)) || mid(a) - mid(b); });
    baits.forEach(function (i) {
      var names = geneNames(D.graph.nodes.g[i]);
      if (!names.length || shown.has(names[0])) return;
      var x = X(mid(i));
      if (placed.some(function (p) { return Math.abs(p - x) < minGap; })) return;
      placed.push(x);
      shown.add(names[0]);
      var t = svg('text', { x: 0, y: 0, 'font-size': 10.5, 'font-style': 'italic', 'text-anchor': 'end',
                            transform: 'translate(' + (x + 3.5).toFixed(1) + ',' + (yC + hC + 6) + ') rotate(-90)',
                            'font-weight': ctx.seeds.indexOf(i) >= 0 ? 700 : 400 });
      t.textContent = names[0].length > 11 ? names[0].slice(0, 10) + '…' : names[0];
      s.appendChild(t);
    });
    // genomic axis in Mb
    var yAx = yC + hC + labelH + 4;
    s.appendChild(svg('line', { x1: L, x2: R, y1: yAx, y2: yAx, class: 'tick' }));
    var span = (x1 - x0) || 1, stepMb = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5].filter(function (st) { return span / 1e6 / st <= (narrow ? 4 : 7); })[0] || 5;
    for (var t = Math.ceil((x0 - pad) / 1e6 / stepMb) * stepMb; t * 1e6 <= x1 + pad; t += stepMb) {
      var xt = X(t * 1e6);
      s.appendChild(svg('line', { x1: xt, x2: xt, y1: yAx, y2: yAx + 4, class: 'tick' }));
      var tt = svg('text', { x: xt, y: yAx + 16, 'font-size': 10.5, 'text-anchor': 'middle', class: 'muted' });
      tt.textContent = t.toFixed(stepMb < 0.1 ? 2 : stepMb < 1 ? 1 : 0);
      s.appendChild(tt);
    }
    var cl = svg('text', { x: (L + R) / 2, y: yAx + 32, 'font-size': 10.5, 'text-anchor': 'middle', class: 'muted' });
    cl.textContent = chrom(ctx.seeds[0]) + ' (Mb)';
    s.appendChild(cl);
    host.appendChild(s);
    host.appendChild(legend([
      ['line', 'var(--shared)', 'interaction of the condition / shared'],
      ['line', 'var(--gained)', 'gained (only in ' + M.names[1] + ')'],
      ['line', 'var(--lost)', 'lost (only in ' + M.names[0] + ')'],
      ['seed', '', 'promoter of the gene'], ['bait', '', 'baited promoter'], ['expl', '', 'carries the explanation']]));
    host.appendChild(el('figcaption', null, '<span class="label">Arc view.</span> Interactions around <x-g>' + esc(ctx.gene) + '</x-g> up to ' + ctx.k +
      (ctx.k === 1 ? ' step' : ' steps') + ' in ' + esc(M.names[0]) + ' (top), in ' + esc(M.names[1]) + ' (middle) and their change (bottom: gained in ' + esc(M.names[1]) +
      ' above the axis in red, lost below in blue, shared in grey). Line width scales with interaction strength; ' + sub.nodes.length + ' nodes and ' + sub.edges.length + ' interactions' +
      (sub.reduced ? ', reduced from ' + int(sub.total) + ' nodes to the strongest partners and the paths to the nodes that carry the explanation' : '') + '.'));
  }

  function edgeTip(e, ctx) {
    var E = D.graph.edges, M = D.meta;
    var st = STATUS[E.s[e]];
    return '<div class="t">' + (st === 'shared' ? 'Shared interaction' : st === 'gained' ? 'Gained in ' + esc(M.names[1]) : 'Lost in ' + esc(M.names[1])) + '</div>' +
      '<div>' + capital(place(E.u[e], ctx.ref)) + '</div><div>' + capital(place(E.v[e], ctx.ref)) + '</div>' +
      '<div class="m">Strength ' + esc(M.names[0]) + ' ' + sig(E.a[e], 2) + ', ' + esc(M.names[1]) + ' ' + sig(E.b[e], 2) + ' (score / threshold); ' +
      distance(Math.abs(mid(E.u[e]) - mid(E.v[e]))) + ' apart</div>';
  }

  function legend(items) {
    var box = el('div', { class: 'legend' });
    items.forEach(function (it) {
      var s = svg('svg', { width: it[0] === 'gradient' ? 120 : 26, height: 14, 'aria-hidden': 'true' });
      if (it[0] === 'line') s.appendChild(svg('line', { x1: 2, x2: 24, y1: 7, y2: 7, stroke: it[1], 'stroke-width': 2.4 }));
      else if (it[0] === 'width') {
        s.appendChild(svg('line', { x1: 2, x2: 11, y1: 7, y2: 7, class: 'e-cond', 'stroke-width': 0.8 }));
        s.appendChild(svg('line', { x1: 14, x2: 24, y1: 7, y2: 7, class: 'e-cond', 'stroke-width': 3.6 }));
      }
      else if (it[0] === 'seed') s.appendChild(svg('rect', { x: 8, y: 2.5, width: 9, height: 9, class: 'seed-mark' }));
      else if (it[0] === 'bait') s.appendChild(svg('rect', { x: 10.5, y: 4.5, width: 5, height: 5, class: 'bait-mark' }));
      else if (it[0] === 'expl') s.appendChild(svg('circle', { cx: 13, cy: 7, r: 4, class: 'expl-mark' }));
      else if (it[0] === 'square') s.appendChild(svg('rect', { x: 8, y: 2.5, width: 9, height: 9, fill: 'var(--empty)', stroke: 'var(--ink)', 'stroke-width': it[1] || 0.8 }));
      else if (it[0] === 'circle') s.appendChild(svg('circle', { cx: 13, cy: 7, r: 4.5, fill: 'var(--empty)', stroke: it[1] || 'var(--ink)', 'stroke-width': it[3] || 0.8 }));
      else if (it[0] === 'gradient') {
        var id = 'grad' + Math.random().toString(36).slice(2, 8), defs = svg('defs'), g = svg('linearGradient', { id: id });
        [-1, -0.5, 0, 0.5, 1].forEach(function (v, k) { g.appendChild(svg('stop', { offset: (k * 25) + '%', 'stop-color': diverging(v * (D.meta.maxLog2 || 2)) })); });
        defs.appendChild(g);
        s.appendChild(defs);
        s.appendChild(svg('rect', { x: 1, y: 2, width: 118, height: 10, fill: 'url(#' + id + ')', stroke: 'var(--rule)' }));
      }
      var span = el('span', { class: 'item' });
      span.appendChild(s);
      span.appendChild(el('span', null, it[2]));
      box.appendChild(span);
    });
    return box;
  }

  // neighbourhood graph ---------------------------------------------------------------------------------------------
  function drawGraph(host, ctx) {
    var M = D.meta, E = D.graph.edges, N = D.graph.nodes, sub = ctx.sub;
    var grid = el('div', { class: 'graph-grid' });
    var canvas = el('div'), panel = el('div', { class: 'node-panel', 'aria-live': 'polite' });
    grid.appendChild(canvas);
    grid.appendChild(panel);
    host.appendChild(grid);
    var W = Math.max(280, Math.min(canvas.clientWidth || 700, 900));
    var H = Math.round(W < 520 ? W * 1.15 : Math.max(300, Math.min(W * 0.78, 620)));
    var pos = layout(sub, W, H);
    var s = svg('svg', { viewBox: '0 0 ' + W + ' ' + H, width: W, height: H, role: 'img',
                         'aria-label': 'Neighbourhood graph of ' + ctx.gene + ' up to ' + ctx.k + ' steps' });
    var edgeLayer = svg('g'), nodeLayer = svg('g'), labelLayer = svg('g');
    s.appendChild(edgeLayer); s.appendChild(nodeLayer); s.appendChild(labelLayer);
    var lines = [];
    sub.edges.forEach(function (e) {
      var st = E.s[e];
      if (state.mode === 'A' && st === 1) return;
      if (state.mode === 'B' && st === 2) return;
      var width = edgeWidth(state.mode === 'change' ? edgeStrength(e) : edgeStrength(e, state.mode));
      var a = pos[E.u[e]], b = pos[E.v[e]];
      var line = svg('line', { x1: a[0].toFixed(1), y1: a[1].toFixed(1), x2: b[0].toFixed(1), y2: b[1].toFixed(1),
                               class: 'gline ' + (state.mode === 'change' ? 'e-' + STATUS[st] : 'e-cond'), 'stroke-width': width.toFixed(2) });
      line._e = e;
      edgeLayer.appendChild(line);
      lines.push(line);
    });
    var marks = {};
    sub.nodes.forEach(function (i) {
      var p = pos[i], bait = !!N.b[i], seed = ctx.seeds.indexOf(i) >= 0;
      var g = svg('g', { class: 'node' + (seed ? ' seed' : '') + (ctx.explSet.has(i) ? ' explaining' : ''), tabindex: 0, role: 'button',
                         'aria-label': 'node ' + coord(i) });
      var fill = bait ? diverging(N.r[i]) : diverging(markChange(i, ctx.mark));
      var shape = bait ? svg('rect', { x: p[0] - (seed ? 6.5 : 5), y: p[1] - (seed ? 6.5 : 5), width: seed ? 13 : 10, height: seed ? 13 : 10, class: 'shape', fill: fill })
                       : svg('circle', { cx: p[0], cy: p[1], r: 4.6, class: 'shape', fill: fill });
      g.appendChild(shape);
      g.addEventListener('mousemove', function (evt) { showTip(evt, nodeSummary(i, ctx)); highlight(i); });
      g.addEventListener('mouseleave', function () { hideTip(); highlight(null); });
      g.addEventListener('click', function () { selectNode(i, ctx); });
      g.addEventListener('keydown', function (evt) { if (evt.key === 'Enter' || evt.key === ' ') { evt.preventDefault(); selectNode(i, ctx); } });
      nodeLayer.appendChild(g);
      marks[i] = g;
    });
    // labels: the gene, promoters that carry the explanation, then the largest expression changes
    var labelled = sub.nodes.filter(function (i) { return N.b[i] && geneNames(N.g[i]).length; });
    labelled.sort(function (a, b) {
      return (ctx.seeds.indexOf(b) >= 0) - (ctx.seeds.indexOf(a) >= 0) || ctx.explSet.has(b) - ctx.explSet.has(a) ||
        (Math.abs(N.r[b] || 0) - Math.abs(N.r[a] || 0));
    });
    var boxes = [], seen = new Set(), maxLabels = W < 500 ? 10 : 22;
    labelled.forEach(function (i) {
      if (boxes.length >= maxLabels) return;
      var name = geneNames(N.g[i])[0];
      if (seen.has(name)) return;
      var p = pos[i], w = name.length * 6.4 + 4, h = 13;
      var candidates = [[p[0] + 8, p[1] - 6], [p[0] - 8 - w, p[1] - 6], [p[0] - w / 2, p[1] - 20], [p[0] - w / 2, p[1] + 8]];
      var spot = candidates.filter(function (c) {
        return c[0] >= 2 && c[0] + w <= W - 2 && c[1] >= 1 && c[1] + h <= H - 1 &&
          !boxes.some(function (b) { return c[0] < b[0] + b[2] && c[0] + w > b[0] && c[1] < b[1] + b[3] && c[1] + h > b[1]; });
      })[0];
      if (!spot) return;
      seen.add(name);
      boxes.push([spot[0], spot[1], w, h]);
      var t = svg('text', { x: spot[0].toFixed(1), y: (spot[1] + 10).toFixed(1), 'font-size': 11, 'font-style': 'italic', class: 'halo',
                            'font-weight': ctx.seeds.indexOf(i) >= 0 ? 700 : 400 });
      t.textContent = name;
      t.style.pointerEvents = 'none';
      labelLayer.appendChild(t);
    });
    var ring = svg('circle', { r: 10, class: 'ring' });
    ring.style.display = 'none';
    s.appendChild(ring);
    canvas.appendChild(s);
    function highlight(i) {
      lines.forEach(function (l) {
        var e = l._e;
        l.classList.toggle('dim', i !== null && E.u[e] !== i && E.v[e] !== i);
      });
    }
    ctx.panel = panel;
    ctx.ring = ring;
    ctx.pos = pos;
    selectNode(state.selected !== null && sub.dist.has(state.selected) ? state.selected : null, ctx, true);
    var counts = [0, 0, 0];
    sub.edges.forEach(function (e) { counts[E.s[e]]++; });
    canvas.appendChild(legend([
      ['square', 0.8, 'baited promoter'], ['circle', 'var(--ink)', 'other end'],
      ['square', 2.6, 'promoter of the gene'], ['circle', 'var(--explaining)', 'carries the explanation', 2.4],
      ['gradient', '', 'log2 ' + M.names[1] + ' / ' + M.names[0] + ', ' + MINUS + sig(M.maxLog2, 2) + ' to +' + sig(M.maxLog2, 2) + ': expression (squares), ' + ctx.mark + ' (circles)']
    ].concat(state.mode === 'change' ? [['line', 'var(--gained)', 'gained'], ['line', 'var(--lost)', 'lost'], ['line', 'var(--shared)', 'shared']] : [['width', '', 'width by strength in ' + (state.mode === 'A' ? M.names[0] : M.names[1])]])));
    var named = [];
    [1, 2, 3].forEach(function (step) {
      var items = ctx.expl.filter(function (x) { return x[1] === step; }).slice(0, 3)
        .map(function (x) { return place(x[0], ctx.ref) + ' (' + esc(x[2]) + ' ' + log2txt(x[3]) + ')'; });
      if (items.length) named.push((step === 1 ? 'direct partners' : step + ' steps away') + ': ' + joinWords(items));
    });
    host.appendChild(el('figcaption', null, '<span class="label">Neighbourhood graph.</span> <x-g>' + esc(ctx.gene) + '</x-g> and the nodes up to ' + ctx.k +
      (ctx.k === 1 ? ' step' : ' steps') + ' away (' + (state.mode === 'change' ? counts[1] + ' gained, ' + counts[2] + ' lost, ' + counts[0] + ' shared interactions' :
        'interactions in ' + esc(state.mode === 'A' ? M.names[0] : M.names[1])) + '). Squares: baited promoters, coloured by the expression change; circles: other ends, coloured by the change of ' +
      esc(ctx.mark) + '. Black outline: the promoter of the gene; purple outline: nodes that carry the explanation' + (named.length ? ' (' + named.join('; ') + ')' : '') +
      '. Line width scales with interaction strength.' + (sub.reduced ? ' Reduced from ' + int(sub.total) + ' nodes to the strongest partners and the paths to the nodes that carry the explanation.' : '')));
    linkGenes(host);
  }

  function selectNode(i, ctx, quiet) {
    state.selected = i;
    if (ctx.panel) nodePanel(ctx.panel, i, ctx);
    if (ctx.ring) {
      if (i !== null && ctx.pos && ctx.pos[i]) {
        ctx.ring.setAttribute('cx', ctx.pos[i][0]);
        ctx.ring.setAttribute('cy', ctx.pos[i][1]);
        ctx.ring.style.display = '';
      } else ctx.ring.style.display = 'none';
    }
    if (!quiet && ctx.panel && window.innerWidth <= 860) ctx.panel.scrollIntoView({ block: 'nearest', behavior: 'smooth' });
  }

  // triggers and switching events of the gene -----------------------------------------------------------------------
  function columnIndex(table) { var out = {}; table.cols.forEach(function (c, k) { out[c] = k; }); return out; }

  function geneTriggers(page, gene) {
    if (!D.triggers) return;
    var c = columnIndex(D.triggers);
    var rows = D.triggers.rows.filter(function (r) { return String(r[c.gene]).split(',').indexOf(gene) >= 0; });
    if (!rows.length) return;
    section(page, 'Candidate triggers', 'triggers');
    page.appendChild(el('p', { class: 'note' }, 'Transcription factors with a motif at the nodes that carry the explanation of <x-g>' + esc(gene) +
      '</x-g> (egaTriggers), the best ' + D.triggers.maxRank + ' by rank.'));
    var host = el('div');
    page.appendChild(host);
    makeTable(host, { rows: rows, columns: triggerColumns(c, false), sort: 'rank', page: 20, noun: ['factor', 'factors'] });
  }

  function geneSwitching(page, gene) {
    if (!D.switching) return;
    var c = columnIndex(D.switching);
    var rows = D.switching.rows.filter(function (r) {
      return [r[c.genes_lost], r[c.genes_gained]].some(function (t) { return String(t || '').split(',').indexOf(gene) >= 0; });
    });
    if (!rows.length) return;
    section(page, 'Enhancer switching', 'switching');
    page.appendChild(el('p', { class: 'note' }, 'Elements that lose their contact with one promoter and gain one with another between ' + esc(D.meta.names[0]) + ' and ' +
      esc(D.meta.names[1]) + ', where <x-g>' + esc(gene) + '</x-g> is one of the two promoters.'));
    var host = el('div');
    page.appendChild(host);
    makeTable(host, { rows: rows, columns: switchColumns(c), sort: 'rank', page: 20, noun: ['event', 'events'] });
  }

  // ------------------------------------------------------------------------------------------------------------------
  // browse: differential genes

  function renderGenes() {
    var v = view(), T = D.table, c = columnIndex(T), M = D.meta;
    v.appendChild(el('h2', { text: 'Differential genes' }));
    v.appendChild(el('p', { class: 'note' }, 'One row per differential promoter fragment. First level: the nearest level that explains the expression change; ' +
      'explaining levels: every level that explains it. Best q: the smallest q of an explaining step away from the promoter; distance control: a level that explains the gene ' +
      'is significant (q &lt; 0.05) and beyond distance-matched mirror loci. Click a row to open the gene.'));
    var f = { text: '', level: '', dir: '', supported: false, cls: '' };
    var controls = controlsRow(v);
    textControl(controls, 'Filter', 'gene name', function (x) { f.text = x; apply(); });
    var levels = Object.keys(M.levelLabels).sort(function (a, b) { return M.levelOrder[a] - M.levelOrder[b]; });
    selectControl(controls, 'Explained at', [['', 'any level']].concat(levels.map(function (lv) { return [lv, levelShort(lv)]; })).concat([['none', 'no level']]), '',
      function (x) { f.level = x; apply(); });
    selectControl(controls, 'Direction', [['', 'up and down'], ['up', 'up in ' + M.names[1]], ['down', 'down in ' + M.names[1]]], '', function (x) { f.dir = x; apply(); });
    var classSet = {};
    T.rows.forEach(function (r) { if (r[c.class_a]) classSet[r[c.class_a]] = 1; if (r[c.class_b]) classSet[r[c.class_b]] = 1; });
    var classList = Object.keys(classSet).sort();
    if (classList.length) selectControl(controls, 'Promoter class', [['', 'any']].concat(classList.map(function (x) { return [x, x]; })), '', function (x) { f.cls = x; apply(); });
    checkControl(controls, 'passes the distance control', function (x) { f.supported = x; apply(); });
    var host = el('div');
    v.appendChild(host);
    var table = makeTable(host, {
      rows: T.rows, sort: 'rna_abs', sortDir: -1, page: 100, noun: ['fragment', 'fragments'],
      onRow: function (r) { setHash({ gene: r[c.gene] }); },
      columns: [
        { id: 'gene', label: 'Gene', value: function (r) { return r[c.gene]; }, render: function (r) {
          var others = r[c.genes].split(',').filter(function (g) { return g && g !== r[c.gene]; });
          return geneTag(r[c.gene]) + (others.length ? ' <span class="note">(with ' + others.slice(0, 3).map(geneTag).join(', ') + (others.length > 3 ? ', ...' : '') + ')</span>' : '');
        } },
        { id: 'rna', label: 'log2 FC', num: true, value: function (r) { return r[c.rna]; }, render: function (r) { return '<span class="' + (r[c.rna] > 0 ? 'dir-up' : 'dir-down') + '">' + log2txt(r[c.rna]) + '</span>'; } },
        { id: 'rna_abs', label: '|log2 FC|', num: true, value: function (r) { return Math.abs(r[c.rna]); }, render: function (r) { return sig(Math.abs(r[c.rna]), 2); } },
        { id: 'first', label: 'First level', value: function (r) { return r[c.first] ? M.levelOrder[r[c.first]] : 999; }, render: function (r) { return r[c.first] ? esc(levelShort(r[c.first])) : '<span class="note">none</span>'; } },
        { id: 'levels', label: 'Explaining levels', value: function (r) { return r[c.levels].length; }, render: function (r) { return esc(r[c.levels].map(levelShort).join(', ')); } },
        { id: 'q', label: 'Best q', num: true, value: function (r) { return r[c.q]; }, render: function (r) { return isnum(r[c.q]) ? pq(r[c.q]) + ' <span class="note">(' + esc(levelShort(r[c.q_level])) + ')</span>' : ''; } },
        { id: 'supported', label: 'Distance control', value: function (r) { return -r[c.supported]; }, render: function (r) { return r[c.supported] ? '<span class="yes">passed</span>' : ''; } },
        { id: 'class', label: 'Promoter class', value: function (r) { return r[c.class_a] + r[c.class_b]; }, render: function (r) {
          return r[c.class_a] === r[c.class_b] ? esc(r[c.class_a]) : esc(r[c.class_a] || 'n.a.') + ' → ' + esc(r[c.class_b] || 'n.a.');
        } }
      ]
    });
    function apply() {
      table.filter(function (r) {
        if (f.text && r[c.genes].toLowerCase().indexOf(f.text) < 0) return false;
        if (f.level === 'none' && r[c.levels].length) return false;
        if (f.level && f.level !== 'none' && r[c.levels].indexOf(f.level) < 0) return false;
        if (f.dir === 'up' && !(r[c.rna] > 0)) return false;
        if (f.dir === 'down' && !(r[c.rna] < 0)) return false;
        if (f.supported && !r[c.supported]) return false;
        if (f.cls && r[c.class_a] !== f.cls && r[c.class_b] !== f.cls) return false;
        return true;
      });
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  // browse: communities

  function renderCommunities() {
    var v = view(), C = D.communities, M = D.meta;
    v.appendChild(el('h2', { text: 'Communities' }));
    if (!C) { v.appendChild(el('p', { class: 'note', text: 'The run has no community tables.' })); return; }
    var c = columnIndex(C);
    v.appendChild(el('p', { class: 'note' }, 'Communities of promoters in ' + esc(M.names[0]) + ' matched to those of ' + esc(M.names[1]) +
      ' (side A), and communities of ' + esc(M.names[1]) + ' without a match (side B, formed). Named by the differential member with the largest expression change. ' +
      'Coherence: the expression of the member genes changes more coherently than in random gene sets of the same size. Click a row for its members.'));
    var detail = el('div');
    var f = { text: '', transition: '', side: '', coherent: false, named: false };
    var controls = controlsRow(v);
    textControl(controls, 'Filter', 'gene name', function (x) { f.text = x; apply(); });
    var transitions = {};
    C.rows.forEach(function (r) { transitions[r[c.transition]] = 1; });
    selectControl(controls, 'Transition', [['', 'any']].concat(Object.keys(transitions).sort().map(function (t) { return [t, t]; })), '', function (x) { f.transition = x; apply(); });
    selectControl(controls, 'Side', [['', 'A and B'], ['A', 'A (' + M.names[0] + ')'], ['B', 'B (' + M.names[1] + ')']], '', function (x) { f.side = x; apply(); });
    checkControl(controls, 'coherent (q < 0.1)', function (x) { f.coherent = x; apply(); });
    checkControl(controls, 'with a differential member', function (x) { f.named = x; apply(); });
    v.appendChild(detail);
    var host = el('div');
    v.appendChild(host);
    // member gene names, for the text filter
    var memberNames = {};
    var mem = C.members;
    mem.node.forEach(function (node, k) {
      var names = geneNames(D.graph.nodes.g[node]);
      if (!names.length) return;
      ['a', 'b'].forEach(function (side) {
        var id = mem[side][k];
        if (id === null) return;
        var key = side + id;
        memberNames[key] = (memberNames[key] || '') + ',' + names.join(',').toLowerCase();
      });
    });
    var table = makeTable(host, {
      rows: C.rows, sort: 'size_a', sortDir: -1, page: 100, noun: ['community', 'communities'],
      onRow: function (r) { showCommunity(detail, r, c); detail.scrollIntoView({ block: 'nearest', behavior: 'smooth' }); },
      columns: [
        { id: 'name', label: 'Named by', value: function (r) { return r[c.name] || '~'; }, render: function (r) { return r[c.name] ? geneTag(r[c.name]) : ''; } },
        { id: 'side', label: 'Side', value: function (r) { return r[c.side]; } },
        { id: 'a', label: 'Community ' + M.names[0], num: true, value: function (r) { return r[c.a]; }, render: function (r) { return r[c.a] === null ? '' : String(r[c.a]); } },
        { id: 'size_a', label: 'Promoters', num: true, value: function (r) { return r[c.size_a]; }, render: function (r) { return int(r[c.size_a]); } },
        { id: 'transition', label: 'Transition', value: function (r) { return r[c.transition]; } },
        { id: 'b', label: 'Community ' + M.names[1], num: true, value: function (r) { return r[c.b]; }, render: function (r) { return r[c.b] === null ? '' : String(r[c.b]); } },
        { id: 'size_b', label: 'Promoters', num: true, value: function (r) { return r[c.size_b]; }, render: function (r) { return int(r[c.size_b]); } },
        { id: 'jaccard', label: 'Shared members', num: true, value: function (r) { return r[c.jaccard]; }, render: function (r) { return isnum(r[c.jaccard]) ? (100 * r[c.jaccard]).toFixed(0) + '%' : ''; } },
        { id: 'n_genes', label: 'Measured genes', num: true, value: function (r) { return r[c.n_genes]; }, render: function (r) { return int(r[c.n_genes]); } },
        { id: 'mean_rna', label: 'Mean log2 FC', num: true, value: function (r) { return r[c.mean_rna]; }, render: function (r) { return log2txt(r[c.mean_rna]); } },
        { id: 'q', label: 'Coherence q', num: true, value: function (r) { return r[c.q]; }, render: function (r) { return pq(r[c.q]); } },
        { id: 'rewired', label: 'Rewired', num: true, value: function (r) { return r[c.rewired]; }, render: function (r) { return isnum(r[c.rewired]) ? (100 * r[c.rewired]).toFixed(0) + '%' : ''; } }
      ]
    });
    function apply() {
      table.filter(function (r) {
        if (f.transition && r[c.transition] !== f.transition) return false;
        if (f.side && r[c.side] !== f.side) return false;
        if (f.coherent && !(r[c.q] < 0.1)) return false;
        if (f.named && !r[c.name]) return false;
        if (f.text) {
          var text = (r[c.a] !== null ? memberNames['a' + r[c.a]] || '' : '') + (r[c.b] !== null ? memberNames['b' + r[c.b]] || '' : '');
          if (text.indexOf(f.text) < 0) return false;
        }
        return true;
      });
    }
  }

  function showCommunity(host, r, c) {
    var C = D.communities, M = D.meta, mem = C.members, N = D.graph.nodes;
    var a = r[c.a], b = r[c.b];
    var inA = new Set(), inB = new Set(), diff = new Set();
    mem.node.forEach(function (node, k) {
      if (a !== null && mem.a[k] === a) inA.add(node);
      if (b !== null && mem.b[k] === b) inB.add(node);
      if (mem.d[k]) diff.add(node);
    });
    var all = Array.from(new Set(Array.from(inA).concat(Array.from(inB))));
    var values = all.filter(function (x) { return diff.has(x) && isnum(N.r[x]); }).map(function (x) { return N.r[x]; });
    var mean = values.length ? values.reduce(function (s, x) { return s + x; }, 0) / values.length : 0;
    var dir = mean >= 0 ? 1 : -1;
    function names(list, withValue, limit) {
      var seen = new Set(), out = [];
      list.sort(function (x, y) { return Math.abs(N.r[y] || 0) - Math.abs(N.r[x] || 0); }).forEach(function (x) {
        var nm = geneNames(N.g[x])[0];
        if (!nm || seen.has(nm)) return;
        seen.add(nm);
        out.push(geneTag(nm) + (withValue && isnum(N.r[x]) ? ' ' + log2txt(N.r[x]) : ''));
      });
      var more = out.length - limit;
      return { text: out.slice(0, limit).join(', ') + (more > 0 ? ' and ' + more + ' more' : ''), n: out.length };
    }
    var same = names(all.filter(function (x) { return diff.has(x) && Math.sign(N.r[x]) === dir; }), true, 30);
    var opposite = names(all.filter(function (x) { return diff.has(x) && Math.sign(N.r[x]) === -dir; }), true, 30);
    var left = names(Array.from(inA).filter(function (x) { return !inB.has(x); }), false, 40);
    var joined = names(Array.from(inB).filter(function (x) { return !inA.has(x); }), false, 40);
    var both = Array.from(inA).filter(function (x) { return inB.has(x); }).length;
    var html = '<div class="card"><h4>' + (r[c.name] ? 'Community of ' + geneTag(r[c.name]) : 'Community') + '</h4><p>' +
      (a !== null ? 'Community ' + a + ' of ' + esc(M.names[0]) + ' (' + int(inA.size) + ' promoters)' : 'No community of ' + esc(M.names[0])) + '; ' +
      (b !== null ? 'community ' + b + ' of ' + esc(M.names[1]) + ' (' + int(inB.size) + ' promoters)' : 'no community of ' + esc(M.names[1])) +
      '; ' + int(both) + ' promoters in both; transition: ' + esc(r[c.transition]) + '.' +
      (isnum(r[c.q]) ? ' Expression of ' + int(r[c.n_genes]) + ' measured member genes: mean log2 ' + log2txt(r[c.mean_rna]) + ', P = ' + pq(r[c.p]) + ', q = ' + pq(r[c.q]) + '.' : '') + '</p>' +
      '<p><strong>Same direction (' + (dir > 0 ? 'up' : 'down') + ', ' + same.n + '):</strong> ' + (same.text || 'none') + '</p>' +
      '<p><strong>Opposite direction (' + opposite.n + '):</strong> ' + (opposite.text || 'none') + '</p>' +
      '<p><strong>Left in ' + esc(M.names[1]) + ' (' + left.n + '):</strong> ' + (left.text || 'none') + '</p>' +
      '<p><strong>Joined in ' + esc(M.names[1]) + ' (' + joined.n + '):</strong> ' + (joined.text || 'none') + '</p></div>';
    host.innerHTML = html;
    linkGenes(host);
  }

  // ------------------------------------------------------------------------------------------------------------------
  // browse: enhancer switching

  function switchColumns(c) {
    var M = D.meta;
    function promoter(kind) {
      return function (r) {
        var genes = String(r[c['genes_' + kind]] || '').split(',').filter(Boolean);
        return genes.slice(0, 3).map(geneTag).join(', ') + (genes.length > 3 ? ', ...' : '') +
          (isnum(r[c['rna_log2FC_' + kind]]) ? ' <span class="note">' + log2txt(r[c['rna_log2FC_' + kind]]) + '</span>' : '');
      };
    }
    var cols = [
      { id: 'rank', label: 'Rank', num: true, value: function (r) { return r[c.rank]; } },
      { id: 'element', label: 'Element', value: function (r) { return r[c.element]; }, render: function (r) {
        var i = r[c.element];
        var type = esc(r[c.element_type] || '');
        var genes = String(r[c.element_genes] || '').split(',').filter(Boolean);
        return type + (genes.length ? ' (' + genes.slice(0, 2).map(geneTag).join(', ') + ')' : '') + (i !== null && i !== undefined ? '<br><span class="note mono">' + esc(coord(i)) + '</span>' : '');
      }, cls: 'long' }
    ];
    if (c.class_a !== undefined) cols.push({ id: 'class', label: 'Class ' + M.names[0] + ' → ' + M.names[1], value: function (r) { return r[c.class_a] + r[c.class_b]; },
      render: function (r) { return r[c.class_a] === r[c.class_b] ? esc(r[c.class_a]) : esc(r[c.class_a]) + ' → ' + esc(r[c.class_b]); } });
    cols = cols.concat([
      { id: 'lost', label: 'Contact lost', value: function (r) { return r[c.genes_lost]; }, render: promoter('lost') },
      { id: 'dlost', label: 'Distance', num: true, value: function (r) { return r[c.distance_lost]; }, render: function (r) { return isnum(r[c.distance_lost]) ? distance(r[c.distance_lost]) : ''; } },
      { id: 'gained', label: 'Contact gained', value: function (r) { return r[c.genes_gained]; }, render: promoter('gained') },
      { id: 'dgained', label: 'Distance', num: true, value: function (r) { return r[c.distance_gained]; }, render: function (r) { return isnum(r[c.distance_gained]) ? distance(r[c.distance_gained]) : ''; } },
      { id: 'concordant', label: 'Concordant', value: function (r) { return -r[c.concordant]; }, render: function (r) { return r[c.concordant] ? '<span class="yes">yes</span>' : 'no'; },
        title: 'the gene that gains the contact goes up and the gene that loses it goes down' },
      { id: 'evidence', label: 'Evidence', num: true, value: function (r) { return r[c.evidence]; }, render: function (r) { return sig(r[c.evidence], 3); } }
    ]);
    return cols;
  }

  function renderSwitching() {
    var v = view(), S = D.switching, M = D.meta;
    v.appendChild(el('h2', { text: 'Enhancer switching' }));
    if (!S) { v.appendChild(el('p', { class: 'note', text: 'No enhancer switching events for this comparison (egaSwitching output not found).' })); return; }
    var c = columnIndex(S);
    v.appendChild(el('p', { class: 'note' }, 'Elements that lose their contact with one promoter and gain a contact with another promoter between ' + esc(M.names[0]) + ' and ' + esc(M.names[1]) +
      ', from ' + esc(S.file) + (S.rows.length < S.total ? ' (the ' + int(S.rows.length) + ' best ranked of ' + int(S.total) + ' events)' : ' (' + int(S.total) + ' events)') +
      '. Concordant: the promoter that gains the contact goes up and the one that loses it goes down. Gene names open the gene.'));
    var f = { text: '', concordant: false, type: '', changed: false };
    var controls = controlsRow(v);
    textControl(controls, 'Filter', 'gene name', function (x) { f.text = x; apply(); });
    var types = {};
    S.rows.forEach(function (r) { types[r[c.element_type]] = 1; });
    selectControl(controls, 'Element type', [['', 'any']].concat(Object.keys(types).sort().map(function (t) { return [t, t]; })), '', function (x) { f.type = x; apply(); });
    checkControl(controls, 'concordant only', function (x) { f.concordant = x; apply(); });
    if (c.element_class_changed !== undefined) checkControl(controls, 'element class changes', function (x) { f.changed = x; apply(); });
    var host = el('div');
    v.appendChild(host);
    var table = makeTable(host, { rows: S.rows, columns: switchColumns(c), sort: 'rank', page: 100, noun: ['event', 'events'] });
    function apply() {
      table.filter(function (r) {
        if (f.concordant && !r[c.concordant]) return false;
        if (f.type && r[c.element_type] !== f.type) return false;
        if (f.changed && !r[c.element_class_changed]) return false;
        if (f.text && (String(r[c.genes_lost]) + ',' + String(r[c.genes_gained]) + ',' + String(r[c.element_genes])).toLowerCase().indexOf(f.text) < 0) return false;
        return true;
      });
    }
  }

  // ------------------------------------------------------------------------------------------------------------------
  // browse: candidate triggers

  function triggerColumns(c, withGene) {
    var cols = [];
    if (withGene) {
      cols.push({ id: 'gene', label: 'Gene', value: function (r) { return r[c.gene]; }, render: function (r) { return String(r[c.gene]).split(',').map(geneTag).join(', '); } });
      cols.push({ id: 'gene_log2FC', label: 'log2 FC', num: true, value: function (r) { return r[c.gene_log2FC]; }, render: function (r) { return log2txt(r[c.gene_log2FC]); } });
      cols.push({ id: 'explained_at', label: 'Explained at', value: function (r) { return r[c.explained_at]; } });
    }
    return cols.concat([
      { id: 'rank', label: 'Rank', num: true, value: function (r) { return r[c.rank]; } },
      { id: 'tf', label: 'Factor', value: function (r) { return r[c.tf]; }, render: function (r) { return D.genes[r[c.tf]] ? geneTag(r[c.tf]) : esc(r[c.tf]); } },
      { id: 'relation', label: 'Relation', value: function (r) { return r[c.relation]; } },
      { id: 'evidence', label: 'Evidence', value: function (r) { return r[c.evidence]; } },
      { id: 'tf_log2FC', label: 'Factor log2 FC', num: true, value: function (r) { return r[c.tf_log2FC]; }, render: function (r) { return log2txt(r[c.tf_log2FC]); } },
      { id: 'motif', label: 'Motif', value: function (r) { return r[c.motif]; } },
      { id: 'steps', label: 'Motif at', value: function (r) { return r[c.steps_with_motif]; }, render: function (r) { return esc(String(r[c.steps_with_motif] || '').replace(/;/g, ', ')); } },
      { id: 'p', label: 'Motif P', num: true, value: function (r) { return r[c.gene_motif_p]; }, render: function (r) { return pq(r[c.gene_motif_p]); } },
      { id: 'score', label: 'Score', num: true, value: function (r) { return r[c.rank_score]; }, render: function (r) { return sig(r[c.rank_score], 3); } }
    ]);
  }

  function renderTriggers() {
    var v = view(), T = D.triggers;
    v.appendChild(el('h2', { text: 'Candidate triggers' }));
    if (!T) { v.appendChild(el('p', { class: 'note', text: 'No candidate triggers for this comparison (egaTriggers output not found).' })); return; }
    var c = columnIndex(T);
    v.appendChild(el('p', { class: 'note' }, 'Transcription factors whose motifs occur at the nodes that carry the explanation of a gene, ranked per gene (egaTriggers, ' + esc(T.file) +
      '); the best ' + T.maxRank + ' per gene, ' + int(T.rows.length) + ' of ' + int(T.total) + ' rows. Relation: whether the change of the factor fits its known role as activator or repressor.'));
    var f = { gene: '', tf: '', relation: '', rank: T.maxRank };
    var controls = controlsRow(v);
    textControl(controls, 'Filter', 'gene name', function (x) { f.gene = x; apply(); });
    textControl(controls, 'Factor', 'factor name', function (x) { f.tf = x; apply(); });
    var relations = {};
    T.rows.forEach(function (r) { relations[r[c.relation]] = 1; });
    selectControl(controls, 'Relation', [['', 'any']].concat(Object.keys(relations).sort().map(function (t) { return [t, t]; })), '', function (x) { f.relation = x; apply(); });
    selectControl(controls, 'Rank up to', [1, 3, 5, 10].filter(function (x) { return x <= T.maxRank; }).map(function (x) { return [x, String(x)]; }).concat([[T.maxRank, String(T.maxRank)]])
      .filter(function (x, k, all) { return all.findIndex(function (y) { return y[0] === x[0]; }) === k; }), T.maxRank, function (x) { f.rank = parseInt(x, 10); apply(); });
    var host = el('div');
    v.appendChild(host);
    var table = makeTable(host, { rows: T.rows, columns: triggerColumns(c, true), sort: 'score', sortDir: -1, page: 100, noun: ['row', 'rows'] });
    function apply() {
      table.filter(function (r) {
        if (f.gene && String(r[c.gene]).toLowerCase().indexOf(f.gene) < 0) return false;
        if (f.tf && String(r[c.tf]).toLowerCase().indexOf(f.tf) < 0) return false;
        if (f.relation && r[c.relation] !== f.relation) return false;
        if (r[c.rank] > f.rank) return false;
        return true;
      });
    }
  }

  // redraw the figures at a new width
  var lastWidth = window.innerWidth, timer = null;
  window.addEventListener('resize', function () {
    if (Math.abs(window.innerWidth - lastWidth) < 40) return;
    lastWidth = window.innerWidth;
    clearTimeout(timer);
    timer = setTimeout(function () { if (D && state.view === 'gene' && document.getElementById('gene-page')) drawFigures(); }, 200);
  });
  if (window.matchMedia) {
    var mq = window.matchMedia('(prefers-color-scheme: dark)');
    var redraw = function () { DIVERGING = null; if (D && state.view === 'gene' && document.getElementById('gene-page')) drawFigures(); };
    if (mq.addEventListener) mq.addEventListener('change', redraw);
  }

  load();
})();
