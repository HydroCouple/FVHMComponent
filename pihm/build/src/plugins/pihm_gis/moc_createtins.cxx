/****************************************************************************
** Meta object code from reading C++ file 'createtins.h'
**
** Created: Tue Jul 27 23:38:15 2010
**      by: The Qt Meta Object Compiler version 59 (Qt 4.3.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../../src/plugins/pihm_gis/DomainDecomposition/CreateTINs/createtins.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'createtins.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 59
#error "This file was generated using the moc from 4.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

static const uint qt_meta_data_createTINsDlg[] = {

 // content:
       1,       // revision
       0,       // classname
       0,    0, // classinfo
       5,   10, // methods
       0,    0, // properties
       0,    0, // enums/sets

 // slots: signature, parameters, type, tag, flags
      15,   14,   14,   14, 0x0a,
      27,   14,   14,   14, 0x0a,
      40,   14,   14,   14, 0x0a,
      52,   14,   14,   14, 0x0a,
      58,   14,   14,   14, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_createTINsDlg[] = {
    "createTINsDlg\0\0eleBrowse()\0nodeBrowse()\0"
    "shpBrowse()\0run()\0help()\0"
};

const QMetaObject createTINsDlg::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_createTINsDlg,
      qt_meta_data_createTINsDlg, 0 }
};

const QMetaObject *createTINsDlg::metaObject() const
{
    return &staticMetaObject;
}

void *createTINsDlg::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_createTINsDlg))
	return static_cast<void*>(const_cast< createTINsDlg*>(this));
    return QDialog::qt_metacast(_clname);
}

int createTINsDlg::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        switch (_id) {
        case 0: eleBrowse(); break;
        case 1: nodeBrowse(); break;
        case 2: shpBrowse(); break;
        case 3: run(); break;
        case 4: help(); break;
        }
        _id -= 5;
    }
    return _id;
}