/****************************************************************************
** Meta object code from reading C++ file 'helpdialog.h'
**
** Created: Tue Feb 20 19:06:08 2007
**      by: The Qt Meta Object Compiler version 59 (Qt 4.1.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "pihmLIBS/helpDialog/helpdialog.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'helpdialog.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.1.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_helpDialog[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets

       0        // eod
};

static const char qt_meta_stringdata_helpDialog[] = {
    "helpDialog\0"
};

const QMetaObject helpDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_helpDialog,
      qt_meta_data_helpDialog, 0 }
};

const QMetaObject *helpDialog::metaObject() const
{
    return &staticMetaObject;
}

void *helpDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_helpDialog))
	return static_cast<void*>(const_cast<helpDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int helpDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    return _id;
}
