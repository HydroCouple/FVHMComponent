/*!
 *  \file    fvhmcompopnentinfo.h
 *  \author  Caleb Amoa Buahin <caleb.buahin@gmail.com>
 *  \version 1.0.0.0
 *  \section Description
 * Main code for the FVHMComponent model.
 *  \section License
 *  fvhmcompopnentinfo.h, associated files and libraries are free software;
 *  you can redistribute it and/or modify it under the terms of the
 *  Lesser GNU General Public License as published by the Free Software Foundation;
 *  either version 3 of the License, or (at your option) any later version.
 *  fvhmcompopnentinfo.h its associated files is distributed in the hope that it will be useful,
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
#include <matrix.h>
#include <vector>
#include <list>


class FVHMCOMPONENT_EXPORT SparseMatrix
{

  public:

    SparseMatrix(int row, int maxColumn);

    ~SparseMatrix();

    const double& operator()(int row, int column) const;

    double& operator()(int row, int column);

    void appendValue(int row, int column, double value);

    void setValue(int row, int column, double value);

    void getDataByRow(double data[]) const;

    void getColumnIndexes(int indexes[]) const;

    int *rows() const { return m_rows  ;}

    int *colsPerRow()  const {return m_ncols;}

    int getColumnCount(int row) const;

    Matrix<double> toMatrix() const;

    int maxColumCount() const;

    int numRows() const ;

    void print() const;

    int getDataSize() const;

  private:
    int m_rowCount = 0, m_maxColCount = 0;
    int m_startColumn = 0, m_startRow = 0;
    int *m_rows;
    int *m_ncols;
    int **m_colIndexes;
    double **m_values;
    double m_zero = 0.0;
    double m_placeH = 0.0;
};



#endif // FVHMMATRIX_H
