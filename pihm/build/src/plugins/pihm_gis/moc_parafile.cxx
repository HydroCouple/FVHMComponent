/****************************************************************************
** Meta object code from reading C++ file 'parafile.h'
**
** Created: Tue Jul 27 23:38:15 2010
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../src/plugins/pihm_gis/DataModelLoader/ParaFile/parafile.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'parafile.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_paraFileDlg[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      13,   12,   12,   12, 0x0a,
      26,   12,   12,   12, 0x0a,
      32,   12,   12,   12, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_paraFileDlg[] = {
    "paraFileDlg\0\0paraBrowse()\0run()\0help()\0"
};

const QMetaObject paraFileDlg::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_paraFileDlg,
      qt_meta_data_paraFileDlg, 0 }
};

const QMetaObject *paraFileDlg::metaObject() const
{
    return &staticMetaObject;
}

void *paraFileDlg::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_paraFileDlg))
	return static_cast<void*>(const_cast< paraFileDlg*>(this));
    return QDialog::qt_metacast(_clname);
}

int paraFileDlg::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: paraBrowse(); break;
        case 1: run(); break;
        case 2: help(); break;
        }
        _id -= 3;
    }
    return _id;
}
