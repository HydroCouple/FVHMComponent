/********************************************************************************
** Form generated from reading ui file 'runnalldomain.ui'
**
** Created: Tue Jul 27 23:38:17 2010
**      by: Qt User Interface Compiler version 4.3.2
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_RUNNALLDOMAIN_H
#define UI_RUNNALLDOMAIN_H

#include <Qt3Support/Q3TextBrowser>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>

class Ui_RunnAllDomain
{
public:
    QGroupBox *groupBox;
    QPushButton *pushButtonConstrain;
    QLineEdit *lineEditConstrain;
    QLabel *label;
    QGroupBox *groupBox_2;
    QPushButton *pushButtonMesh;
    QLineEdit *lineEditMesh;
    QLabel *label_2;
    QGroupBox *groupBox_3;
    QCheckBox *checkBoxAngle;
    QLineEdit *lineEditAngle;
    QCheckBox *checkBoxArea;
    QLineEdit *lineEditArea;
    QCheckBox *checkBoxOthers;
    QLineEdit *lineEditOthers;
    QGroupBox *groupBox_4;
    QPushButton *pushButtonRun;
    QPushButton *pushButtonClose;
    QPushButton *pushButtonHelp;
    Q3TextBrowser *textBrowser;

    void setupUi(QDialog *RunnAllDomain)
    {
    if (RunnAllDomain->objectName().isEmpty())
        RunnAllDomain->setObjectName(QString::fromUtf8("RunnAllDomain"));
    RunnAllDomain->resize(600, 470);
    groupBox = new QGroupBox(RunnAllDomain);
    groupBox->setObjectName(QString::fromUtf8("groupBox"));
    groupBox->setGeometry(QRect(15, 5, 571, 80));
    pushButtonConstrain = new QPushButton(groupBox);
    pushButtonConstrain->setObjectName(QString::fromUtf8("pushButtonConstrain"));
    pushButtonConstrain->setGeometry(QRect(460, 36, 98, 32));
    lineEditConstrain = new QLineEdit(groupBox);
    lineEditConstrain->setObjectName(QString::fromUtf8("lineEditConstrain"));
    lineEditConstrain->setGeometry(QRect(150, 40, 281, 22));
    label = new QLabel(groupBox);
    label->setObjectName(QString::fromUtf8("label"));
    label->setGeometry(QRect(14, 40, 121, 17));
    groupBox_2 = new QGroupBox(RunnAllDomain);
    groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
    groupBox_2->setGeometry(QRect(15, 173, 571, 80));
    pushButtonMesh = new QPushButton(groupBox_2);
    pushButtonMesh->setObjectName(QString::fromUtf8("pushButtonMesh"));
    pushButtonMesh->setGeometry(QRect(460, 37, 98, 32));
    lineEditMesh = new QLineEdit(groupBox_2);
    lineEditMesh->setObjectName(QString::fromUtf8("lineEditMesh"));
    lineEditMesh->setGeometry(QRect(150, 40, 281, 22));
    label_2 = new QLabel(groupBox_2);
    label_2->setObjectName(QString::fromUtf8("label_2"));
    label_2->setGeometry(QRect(15, 40, 121, 17));
    groupBox_3 = new QGroupBox(RunnAllDomain);
    groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
    groupBox_3->setGeometry(QRect(15, 90, 571, 80));
    checkBoxAngle = new QCheckBox(groupBox_3);
    checkBoxAngle->setObjectName(QString::fromUtf8("checkBoxAngle"));
    checkBoxAngle->setGeometry(QRect(12, 40, 87, 21));
    lineEditAngle = new QLineEdit(groupBox_3);
    lineEditAngle->setObjectName(QString::fromUtf8("lineEditAngle"));
    lineEditAngle->setEnabled(true);
    lineEditAngle->setGeometry(QRect(80, 40, 61, 22));
    lineEditAngle->setFrame(true);
    checkBoxArea = new QCheckBox(groupBox_3);
    checkBoxArea->setObjectName(QString::fromUtf8("checkBoxArea"));
    checkBoxArea->setGeometry(QRect(200, 36, 91, 31));
    QFont font;
    font.setBold(false);
    font.setWeight(50);
    checkBoxArea->setFont(font);
    lineEditArea = new QLineEdit(groupBox_3);
    lineEditArea->setObjectName(QString::fromUtf8("lineEditArea"));
    lineEditArea->setGeometry(QRect(300, 40, 61, 22));
    checkBoxOthers = new QCheckBox(groupBox_3);
    checkBoxOthers->setObjectName(QString::fromUtf8("checkBoxOthers"));
    checkBoxOthers->setGeometry(QRect(412, 40, 71, 21));
    lineEditOthers = new QLineEdit(groupBox_3);
    lineEditOthers->setObjectName(QString::fromUtf8("lineEditOthers"));
    lineEditOthers->setGeometry(QRect(490, 40, 61, 22));
    groupBox_4 = new QGroupBox(RunnAllDomain);
    groupBox_4->setObjectName(QString::fromUtf8("groupBox_4"));
    groupBox_4->setGeometry(QRect(15, 268, 571, 53));
    pushButtonRun = new QPushButton(groupBox_4);
    pushButtonRun->setObjectName(QString::fromUtf8("pushButtonRun"));
    pushButtonRun->setGeometry(QRect(461, 12, 98, 32));
    pushButtonClose = new QPushButton(groupBox_4);
    pushButtonClose->setObjectName(QString::fromUtf8("pushButtonClose"));
    pushButtonClose->setGeometry(QRect(340, 12, 98, 32));
    pushButtonHelp = new QPushButton(groupBox_4);
    pushButtonHelp->setObjectName(QString::fromUtf8("pushButtonHelp"));
    pushButtonHelp->setGeometry(QRect(10, 12, 98, 32));
    textBrowser = new Q3TextBrowser(RunnAllDomain);
    textBrowser->setObjectName(QString::fromUtf8("textBrowser"));
    textBrowser->setGeometry(QRect(15, 337, 571, 121));

    retranslateUi(RunnAllDomain);
    QObject::connect(checkBoxAngle, SIGNAL(toggled(bool)), lineEditAngle, SLOT(setShown(bool)));
    QObject::connect(checkBoxArea, SIGNAL(toggled(bool)), lineEditArea, SLOT(setShown(bool)));
    QObject::connect(checkBoxOthers, SIGNAL(toggled(bool)), lineEditOthers, SLOT(setShown(bool)));

    QMetaObject::connectSlotsByName(RunnAllDomain);
    } // setupUi

    void retranslateUi(QDialog *RunnAllDomain)
    {
    RunnAllDomain->setWindowTitle(QApplication::translate("RunnAllDomain", "RunnAllDomain", 0, QApplication::UnicodeUTF8));
    groupBox->setTitle(QApplication::translate("RunnAllDomain", "Input", 0, QApplication::UnicodeUTF8));
    pushButtonConstrain->setText(QApplication::translate("RunnAllDomain", "Browse", 0, QApplication::UnicodeUTF8));
    label->setToolTip(QApplication::translate("RunnAllDomain", "Input Constraining layer for Domain Decomposition. Generally, an output shapefile produced after Vector Merge.", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("RunnAllDomain", "Constrain Layer", 0, QApplication::UnicodeUTF8));
    groupBox_2->setTitle(QApplication::translate("RunnAllDomain", "Output", 0, QApplication::UnicodeUTF8));
    pushButtonMesh->setText(QApplication::translate("RunnAllDomain", "Browse", 0, QApplication::UnicodeUTF8));
    label_2->setToolTip(QApplication::translate("RunnAllDomain", "Specify a Shape file where TINs generated by TRIANGLE can be save.", 0, QApplication::UnicodeUTF8));
    label_2->setText(QApplication::translate("RunnAllDomain", "Mesh (TIN) File", 0, QApplication::UnicodeUTF8));
    groupBox_3->setTitle(QApplication::translate("RunnAllDomain", "Options", 0, QApplication::UnicodeUTF8));
    checkBoxAngle->setToolTip(QApplication::translate("RunnAllDomain", "Quality 'Angle' constraint for TIN generation", 0, QApplication::UnicodeUTF8));
    checkBoxAngle->setText(QApplication::translate("RunnAllDomain", "Angle", 0, QApplication::UnicodeUTF8));
    checkBoxArea->setToolTip(QApplication::translate("RunnAllDomain", "Maximum Area Quality constraint for TIN generation.", 0, QApplication::UnicodeUTF8));
    checkBoxArea->setText(QApplication::translate("RunnAllDomain", "Area (m^2)", 0, QApplication::UnicodeUTF8));
    checkBoxOthers->setToolTip(QApplication::translate("RunnAllDomain", "Other TRINAGLE options", 0, QApplication::UnicodeUTF8));
    checkBoxOthers->setText(QApplication::translate("RunnAllDomain", "Others", 0, QApplication::UnicodeUTF8));
    groupBox_4->setTitle(QString());
    pushButtonRun->setText(QApplication::translate("RunnAllDomain", "Run", 0, QApplication::UnicodeUTF8));
    pushButtonClose->setText(QApplication::translate("RunnAllDomain", "Close", 0, QApplication::UnicodeUTF8));
    pushButtonHelp->setText(QApplication::translate("RunnAllDomain", "Help", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(RunnAllDomain);
    } // retranslateUi

};

namespace Ui {
    class RunnAllDomain: public Ui_RunnAllDomain {};
} // namespace Ui

#endif // UI_RUNNALLDOMAIN_H
