# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'image_translation.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(912, 635)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.frame = QtWidgets.QFrame(Dialog)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.splitter = QtWidgets.QSplitter(self.frame)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.leftGraphWidget = ImageView(self.splitter)
        self.leftGraphWidget.setObjectName("leftGraphWidget")
        self.gridLayoutWidget = QtWidgets.QWidget(self.splitter)
        self.gridLayoutWidget.setObjectName("gridLayoutWidget")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.gridLayoutWidget)
        self.gridLayout_2.setContentsMargins(0, 0, 0, 0)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.leftButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.leftButton.sizePolicy().hasHeightForWidth())
        self.leftButton.setSizePolicy(sizePolicy)
        self.leftButton.setObjectName("leftButton")
        self.gridLayout_2.addWidget(self.leftButton, 1, 0, 1, 1)
        self.rightButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.rightButton.sizePolicy().hasHeightForWidth())
        self.rightButton.setSizePolicy(sizePolicy)
        self.rightButton.setObjectName("rightButton")
        self.gridLayout_2.addWidget(self.rightButton, 1, 2, 1, 1)
        self.clockwiseButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.clockwiseButton.setObjectName("clockwiseButton")
        self.gridLayout_2.addWidget(self.clockwiseButton, 4, 2, 1, 1)
        self.downButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.downButton.setObjectName("downButton")
        self.gridLayout_2.addWidget(self.downButton, 2, 1, 1, 1)
        self.counterClockwiseButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.counterClockwiseButton.setObjectName("counterClockwiseButton")
        self.gridLayout_2.addWidget(self.counterClockwiseButton, 4, 0, 1, 1)
        self.upButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(1)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.upButton.sizePolicy().hasHeightForWidth())
        self.upButton.setSizePolicy(sizePolicy)
        self.upButton.setObjectName("upButton")
        self.gridLayout_2.addWidget(self.upButton, 0, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_3.setObjectName("label_3")
        self.gridLayout_2.addWidget(self.label_3, 4, 3, 1, 1)
        self.bandsRefEnergy = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.bandsRefEnergy.setObjectName("bandsRefEnergy")
        self.gridLayout_2.addWidget(self.bandsRefEnergy, 4, 4, 1, 1)
        self.deltaXBox = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.deltaXBox.setObjectName("deltaXBox")
        self.gridLayout_2.addWidget(self.deltaXBox, 0, 4, 1, 1)
        self.label = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 3, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.gridLayoutWidget)
        self.label_2.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout_2.addWidget(self.label_2, 1, 3, 1, 1)
        self.deltaDegBox = QtWidgets.QLineEdit(self.gridLayoutWidget)
        self.deltaDegBox.setObjectName("deltaDegBox")
        self.gridLayout_2.addWidget(self.deltaDegBox, 1, 4, 1, 1)
        self.addBandsButton = QtWidgets.QPushButton(self.gridLayoutWidget)
        self.addBandsButton.setObjectName("addBandsButton")
        self.gridLayout_2.addWidget(self.addBandsButton, 2, 4, 1, 1)
        self.horizontalLayout_2.addWidget(self.splitter)
        self.verticalLayout.addWidget(self.frame)
        self.toolButton = QtWidgets.QToolButton(Dialog)
        self.toolButton.setObjectName("toolButton")
        self.verticalLayout.addWidget(self.toolButton)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.leftButton.setText(_translate("Dialog", "???"))
        self.rightButton.setText(_translate("Dialog", "???"))
        self.clockwiseButton.setText(_translate("Dialog", "???"))
        self.downButton.setText(_translate("Dialog", "???"))
        self.counterClockwiseButton.setText(_translate("Dialog", "???"))
        self.upButton.setText(_translate("Dialog", "???"))
        self.label_3.setText(_translate("Dialog", "Ref Energy"))
        self.label.setText(_translate("Dialog", "??x/??y"))
        self.label_2.setText(_translate("Dialog", "??deg"))
        self.addBandsButton.setText(_translate("Dialog", "Add QE Bands"))
        self.toolButton.setText(_translate("Dialog", "..."))
from pyqtgraph import ImageView
