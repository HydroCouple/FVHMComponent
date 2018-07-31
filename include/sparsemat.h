/*!
 *  \file    sparsemat.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0
 *  \section Description
 * Main code for the FVHMComponent model.
 *  \section License
 *  sparsemat.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU Lesser General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  sparsemat.h its associated files is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.(see <http://www.gnu.org/licenses/> for details)
 *  \date 2014-2016
 *  \pre
 *  \bug
 *  \todo
 *  \warning
 */


#ifndef FVHMMATRIX_H
#define FVHMMATRIX_H

#include "fvhmcomponent_global.h"
#include "matrix.h"
#include <vector>
#include <list>


class FVHMCOMPONENT_EXPORT SparseMatrix
{

  public:

    SparseMatrix(int ilower, int iupper, int columnCount, int maxColSize);

    ~SparseMatrix();

    int ilower() const;

    int iupper() const;

    int rowCount() const ;

    int columnCount() const;

    int maxColumnSize() const;

    const double& operator()(int row, int column) const;

    double& operator()(int row, int column);

    void appendValue(int row, int column, double value);

    void setValue(int row, int column, double value);

    int getDataSize(int ilower, int iupper) const;

    void getValuesByRows(double data[], int ilower, int iupper) const;

    void getColumnIndexes(int indexes[], int ilower, int iupper) const;

    void getRowIndexes(int indexes[], int ilower, int iupper) const;

    void getColsPerRow(int indexes[], int ilower, int iupper)  const;

    int getColumnCount(int row) const;

    Matrix<double> toMatrix() const;

    void serializeRows(int ilower, int iupper, double values[], int &counter) const;

    static void deserialize(const double values[], SparseMatrix* &sparse, double * &b, double *&x, int &counter);

    void print() const;

  private:

    int m_rowCount = 0,
    m_columnCount = 0,
    m_maxColSize = 0,
    m_ilower = 0,
    m_iupper = 0;
//    int *m_rows;
    int *m_ncols;
    int **m_colIndexes;
    double **m_values;
    double m_zero = 0.0;
    double m_placeH = 0.0;
};



#endif // FVHMMATRIX_H
