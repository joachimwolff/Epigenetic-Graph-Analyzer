import sys
from PyQt5 import QtWidgets
import pyqtgraph as pg
import numpy as np
from PyQt5.QtWidgets import QComboBox

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)

        # Create a pop-up dialog to load data
        self.loadDataDialog = QtWidgets.QFileDialog()
        self.loadButton = QtWidgets.QPushButton("Load Data")
        self.loadButton.clicked.connect(self.load_data)

        # Create some widgets
        self.graphWidget1 = pg.PlotWidget()
        self.graphWidget2 = pg.PlotWidget()
        self.dropdown = QComboBox()
        self.dropdown.addItems(['Graph 1', 'Graph 2'])

        # Set a fixed height for the graph widgets
        self.graphWidget1.setFixedHeight(400)
        self.graphWidget2.setFixedHeight(400)

        # Create a layout
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(self.loadButton)
        layout.addWidget(self.dropdown)
        layout.addWidget(self.graphWidget1)
        layout.addWidget(self.graphWidget2)

        # Create a container
        container = QtWidgets.QWidget()
        container.setLayout(layout)

        # Set the container as the central widget
        self.setCentralWidget(container)

        # Generate some data for the graphs
        self.x = None
        self.y = None

        # Connect the dropdown's currentIndexChanged signal to the plot_graph method
        self.dropdown.currentIndexChanged.connect(self.plot_graph)

    def load_data(self):
        # Open the file dialog to select multiple data files
        filepaths, _ = self.loadDataDialog.getOpenFileNames(self, "Load Data", "", "Data Files (*.txt *.csv)")

        # Load the data from the selected files
        if filepaths:
            # Assuming the data is loaded into self.x and self.y
            self.x = np.random.normal(size=1000)
            self.y = np.random.normal(size=(3, 1000))

            # Plot the initial data
            self.plot_graph(0)

    def plot_graph(self, index):
        if self.x is not None and self.y is not None:
            if index == 0:
                self.graphWidget1.clear()
                self.graphWidget1.plot(self.x, self.y[0], pen='r')
            elif index == 1:
                self.graphWidget2.clear()
                self.graphWidget2.plot(self.x, self.y[1], pen='g')

def main():
    app = QtWidgets.QApplication(sys.argv)

    main = MainWindow()
    main.show()

    sys.exit(app.exec_())

if __name__ == '__main__':
    main()