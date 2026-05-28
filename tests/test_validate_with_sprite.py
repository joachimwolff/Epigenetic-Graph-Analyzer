import gzip
import os
import tempfile
import unittest

import networkx as nx

from epigeneticgraphanalyzer.egaValidateWithSPRITE import (
    build_bin_offsets,
    build_graph_triplets,
    build_node_bin_mapping,
    count_target_triplets,
    count_sprite_triplets,
    validate_graph_triplets,
)


class TestValidateWithSPRITE(unittest.TestCase):
    def test_indirect_triplet_supported_by_sprite(self):
        graph = nx.Graph()
        graph.add_node("A", position="chr1:10-30")
        graph.add_node("B", position="chr1:1010-1030")
        graph.add_node("C", position="chr1:2010-2030")
        graph.add_edge("A", "B")
        graph.add_edge("B", "C")

        chrom_sizes = {"chr1": 100000}
        bin_size = 1000
        bin_offsets = build_bin_offsets(chrom_sizes, bin_size)
        node_bin_mapping = build_node_bin_mapping(graph, chrom_sizes, bin_offsets, bin_size)

        with tempfile.TemporaryDirectory() as tmpdir:
            sprite_path = os.path.join(tmpdir, "sprite.tsv.gz")
            with gzip.open(sprite_path, "wt") as handle:
                handle.write("cluster_1\tchr1:10\tchr1:1010\tchr1:2010\n")

            sprite_triplet_counts, _ = count_sprite_triplets(
                sprite_path,
                chrom_sizes,
                bin_offsets,
                bin_size,
            )

        graph_triplets = build_graph_triplets(graph, node_bin_mapping)
        rows, supported = validate_graph_triplets(graph_triplets, sprite_triplet_counts, 1)

        self.assertEqual(len(graph_triplets), 1)
        self.assertEqual(len(rows), 1)
        self.assertEqual(supported, 1)
        self.assertEqual(rows[0]["sprite_triplet_count"], 1)
        self.assertEqual(rows[0]["supported"], 1)

    def test_targeted_mode_counts_only_requested_triplets(self):
        chrom_sizes = {"chr1": 100000}
        bin_size = 1000
        bin_offsets = build_bin_offsets(chrom_sizes, bin_size)
        target_triplets = {(0, 1, 2)}

        with tempfile.TemporaryDirectory() as tmpdir:
            sprite_path = os.path.join(tmpdir, "sprite.tsv.gz")
            with gzip.open(sprite_path, "wt") as handle:
                handle.write("cluster_1\tchr1:10\tchr1:1010\tchr1:2010\tchr1:3010\n")

            counts = count_target_triplets(
                sprite_path,
                chrom_sizes,
                bin_offsets,
                bin_size,
                target_triplets,
            )

        self.assertEqual(counts[(0, 1, 2)], 1)

    def test_targeted_mode_respects_max_cluster_bins(self):
        chrom_sizes = {"chr1": 100000}
        bin_size = 1000
        bin_offsets = build_bin_offsets(chrom_sizes, bin_size)
        target_triplets = {(0, 1, 2)}

        with tempfile.TemporaryDirectory() as tmpdir:
            sprite_path = os.path.join(tmpdir, "sprite.tsv.gz")
            with gzip.open(sprite_path, "wt") as handle:
                handle.write("cluster_1\tchr1:10\tchr1:1010\tchr1:2010\tchr1:3010\n")

            counts = count_target_triplets(
                sprite_path,
                chrom_sizes,
                bin_offsets,
                bin_size,
                target_triplets,
                p_max_cluster_bins=3,
            )

        self.assertEqual(counts[(0, 1, 2)], 0)


if __name__ == "__main__":
    unittest.main()
