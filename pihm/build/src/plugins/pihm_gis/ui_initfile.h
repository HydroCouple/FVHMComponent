/********************************************************************************
** Form generated from reading ui file 'initfile.ui'
**
** Created: Tue Jul 27 23:38:18 2010
**      by: Qt User Interface Compiler version 4.3.2
**
** WARNING! All changes made in this file will be lost when recompiling ui file!
********************************************************************************/

#ifndef UI_INITFILE_H
#define UI_INITFILE_H

#include <Qt3Support/Q3TextBrowser>
#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QDialog>
#include <QtGui/QFrame>
#include <QtGui/QGroupBox>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QPushButton>

class Ui_InitFile
{
public:
    QGroupBox *groupBox;
    QLabel *label;
    QLabel *label_2;
    QLineEdit *lineEditMeshFile;
    QLineEdit *lineEditRivFile;
    QPushButton *pushButton_MeshFile;
    QPushButton *pushButton_RivFile;
    QGroupBox *groupBox_2;
    QLabel *label_3;
    QLabel *label_4;
    QLabel *label_5;
    QLabel *label_6;
    QLabel *label_7;
    QLineEdit *interception;
    QLineEdit *snow;
    QLineEdit *surface;
    QLineEdit *unsaturated;
    QLineEdit *saturated;
    QLabel *label_9;
    QLabel *label_10;
    QLineEdit *river;
    QLineEdit *riverBed;
    QFrame *line;
    QGroupBox *groupBox_3;
    QLabel *label_8;
    QLineEdit *lineEditInitFile;
    QPushButton *pushButton_InitFile;
    QPushButton *pushButton_Help;
    QPushButton *pushButton_Close;
    QPushButton *pushButton_Run;
    Q3TextBrowser *textBrowser;

    void setupUi(QDialog *InitFile)
    {
    if (InitFile->objectName().isEmpty())
        InitFile->setObjectName(QString::fromUtf8("InitFile"));
    InitFile->resize(600, 450);
    groupBox = new QGroupBox(InitFile);
    groupBox->setObjectName(QString::fromUtf8("groupBox"));
    groupBox->setGeometry(QRect(20, 0, 561, 100));
    label = new QLabel(groupBox);
    label->setObjectName(QString::fromUtf8("label"));
    label->setGeometry(QRect(20, 35, 61, 17));
    label_2 = new QLabel(groupBox);
    label_2->setObjectName(QString::fromUtf8("label_2"));
    label_2->setGeometry(QRect(20, 71, 61, 17));
    lineEditMeshFile = new QLineEdit(groupBox);
    lineEditMeshFile->setObjectName(QString::fromUtf8("lineEditMeshFile"));
    lineEditMeshFile->setGeometry(QRect(90, 32, 331, 22));
    lineEditRivFile = new QLineEdit(groupBox);
    lineEditRivFile->setObjectName(QString::fromUtf8("lineEditRivFile"));
    lineEditRivFile->setGeometry(QRect(90, 68, 331, 22));
    pushButton_MeshFile = new QPushButton(groupBox);
    pushButton_MeshFile->setObjectName(QString::fromUtf8("pushButton_MeshFile"));
    pushButton_MeshFile->setGeometry(QRect(430, 28, 113, 32));
    pushButton_RivFile = new QPushButton(groupBox);
    pushButton_RivFile->setObjectName(QString::fromUtf8("pushButton_RivFile"));
    pushButton_RivFile->setGeometry(QRect(431, 64, 113, 32));
    groupBox_2 = new QGroupBox(InitFile);
    groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
    groupBox_2->setGeometry(QRect(20, 101, 561, 121));
    label_3 = new QLabel(groupBox_2);
    label_3->setObjectName(QString::fromUtf8("label_3"));
    label_3->setGeometry(QRect(20, 30, 81, 17));
    label_4 = new QLabel(groupBox_2);
    label_4->setObjectName(QString::fromUtf8("label_4"));
    label_4->setGeometry(QRect(20, 55, 81, 17));
    label_5 = new QLabel(groupBox_2);
    label_5->setObjectName(QString::fromUtf8("label_5"));
    label_5->setGeometry(QRect(220, 40, 61, 17));
    label_6 = new QLabel(groupBox_2);
    label_6->setObjectName(QString::fromUtf8("label_6"));
    label_6->setGeometry(QRect(390, 30, 81, 17));
    label_7 = new QLabel(groupBox_2);
    label_7->setObjectName(QString::fromUtf8("label_7"));
    label_7->setGeometry(QRect(390, 55, 61, 17));
    interception = new QLineEdit(groupBox_2);
    interception->setObjectName(QString::fromUtf8("interception"));
    interception->setGeometry(QRect(110, 26, 61, 22));
    snow = new QLineEdit(groupBox_2);
    snow->setObjectName(QString::fromUtf8("snow"));
    snow->setGeometry(QRect(110, 52, 61, 22));
    surface = new QLineEdit(groupBox_2);
    surface->setObjectName(QString::fromUtf8("surface"));
    surface->setGeometry(QRect(280, 38, 61, 22));
    unsaturated = new QLineEdit(groupBox_2);
    unsaturated->setObjectName(QString::fromUtf8("unsaturated"));
    unsaturated->setGeometry(QRect(480, 26, 61, 22));
    saturated = new QLineEdit(groupBox_2);
    saturated->setObjectName(QString::fromUtf8("saturated"));
    saturated->setGeometry(QRect(480, 52, 61, 22));
    label_9 = new QLabel(groupBox_2);
    label_9->setObjectName(QString::fromUtf8("label_9"));
    label_9->setGeometry(QRect(90, 95, 61, 17));
    label_10 = new QLabel(groupBox_2);
    label_10->setObjectName(QString::fromUtf8("label_10"));
    label_10->setGeometry(QRect(310, 95, 61, 17));
    river = new QLineEdit(groupBox_2);
    river->setObjectName(QString::fromUtf8("river"));
    river->setGeometry(QRect(180, 92, 61, 22));
    riverBed = new QLineEdit(groupBox_2);
    riverBed->setObjectName(QString::fromUtf8("riverBed"));
    riverBed->setGeometry(QRect(410, 92, 61, 22));
    line = new QFrame(groupBox_2);
    line->setObjectName(QString::fromUtf8("line"));
    line->setGeometry(QRect(10, 75, 541, 16));
    line->setLineWidth(0);
    line->setFrameShape(QFrame::HLine);
    line->setFrameShadow(QFrame::Sunken);
    groupBox_3 = new QGroupBox(InitFile);
    groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
    groupBox_3->setGeometry(QRect(20, 219, 561, 71));
    label_8 = new QLabel(groupBox_3);
    label_8->setObjectName(QString::fromUtf8("label_8"));
    label_8->setGeometry(QRect(20, 34, 61, 17));
    lineEditInitFile = new QLineEdit(groupBox_3);
    lineEditInitFile->setObjectName(QString::fromUtf8("lineEditInitFile"));
    lineEditInitFile->setGeometry(QRect(90, 33, 331, 22));
    pushButton_InitFile = new QPushButton(groupBox_3);
    pushButton_InitFile->setObjectName(QString::fromUtf8("pushButton_InitFile"));
    pushButton_InitFile->setGeometry(QRect(432, 30, 113, 32));
    pushButton_Help = new QPushButton(InitFile);
    pushButton_Help->setObjectName(QString::fromUtf8("pushButton_Help"));
    pushButton_Help->setGeometry(QRect(15, 298, 113, 32));
    pushButton_Close = new QPushButton(InitFile);
    pushButton_Close->setObjectName(QString::fromUtf8("pushButton_Close"));
    pushButton_Close->setGeometry(QRect(359, 298, 113, 32));
    pushButton_Run = new QPushButton(InitFile);
    pushButton_Run->setObjectName(QString::fromUtf8("pushButton_Run"));
    pushButton_Run->setGeometry(QRect(475, 298, 113, 32));
    pushButton_Run->setDefault(true);
    textBrowser = new Q3TextBrowser(InitFile);
    textBrowser->setObjectName(QString::fromUtf8("textBrowser"));
    textBrowser->setGeometry(QRect(20, 337, 561, 101));

    retranslateUi(InitFile);

    QMetaObject::connectSlotsByName(InitFile);
    } // setupUi

    void retranslateUi(QDialog *InitFile)
    {
    InitFile->setWindowTitle(QApplication::translate("InitFile", "InitFile", 0, QApplication::UnicodeUTF8));
    groupBox->setTitle(QApplication::translate("InitFile", "Input Files", 0, QApplication::UnicodeUTF8));
    label->setText(QApplication::translate("InitFile", "mesh File", 0, QApplication::UnicodeUTF8));
    label_2->setText(QApplication::translate("InitFile", "riv File", 0, QApplication::UnicodeUTF8));
    pushButton_MeshFile->setText(QApplication::translate("InitFile", "Browse...", 0, QApplication::UnicodeUTF8));
    pushButton_RivFile->setText(QApplication::translate("InitFile", "Browse...", 0, QApplication::UnicodeUTF8));
    groupBox_2->setTitle(QApplication::translate("InitFile", "Uniform Initial Condition (Head in meters) ", 0, QApplication::UnicodeUTF8));
    label_3->setText(QApplication::translate("InitFile", "Interception", 0, QApplication::UnicodeUTF8));
    label_4->setText(QApplication::translate("InitFile", "Snow", 0, QApplication::UnicodeUTF8));
    label_5->setText(QApplication::translate("InitFile", "Surface", 0, QApplication::UnicodeUTF8));
    label_6->setText(QApplication::translate("InitFile", "Unsaturated", 0, QApplication::UnicodeUTF8));
    label_7->setText(QApplication::translate("InitFile", "Saturated", 0, QApplication::UnicodeUTF8));
    interception->setText(QApplication::translate("InitFile", "0", 0, QApplication::UnicodeUTF8));
    snow->setText(QApplication::translate("InitFile", "0", 0, QApplication::UnicodeUTF8));
    surface->setText(QApplication::translate("InitFile", "0", 0, QApplication::UnicodeUTF8));
    unsaturated->setText(QApplication::translate("InitFile", "0", 0, QApplication::UnicodeUTF8));
    saturated->setText(QApplication::translate("InitFile", "4.5", 0, QApplication::UnicodeUTF8));
    label_9->setText(QApplication::translate("InitFile", "River", 0, QApplication::UnicodeUTF8));
    label_10->setText(QApplication::translate("InitFile", "River Bed", 0, QApplication::UnicodeUTF8));
    river->setText(QApplication::translate("InitFile", "0", 0, QApplication::UnicodeUTF8));
    riverBed->setText(QApplication::translate("InitFile", "4.5", 0, QApplication::UnicodeUTF8));
    groupBox_3->setTitle(QApplication::translate("InitFile", "Output", 0, QApplication::UnicodeUTF8));
    label_8->setText(QApplication::translate("InitFile", "init File", 0, QApplication::UnicodeUTF8));
    pushButton_InitFile->setText(QApplication::translate("InitFile", "Browse...", 0, QApplication::UnicodeUTF8));
    pushButton_Help->setText(QApplication::translate("InitFile", "Help", 0, QApplication::UnicodeUTF8));
    pushButton_Close->setText(QApplication::translate("InitFile", "Close", 0, QApplication::UnicodeUTF8));
    pushButton_Run->setText(QApplication::translate("InitFile", "Run", 0, QApplication::UnicodeUTF8));
    Q_UNUSED(InitFile);
    } // retranslateUi

};

namespace Ui {
    class InitFile: public Ui_InitFile {};
} // namespace Ui

#endif // UI_INITFILE_H