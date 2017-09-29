#include "stdafx.h"
#include "sparsemat.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <QString>
#include <QDebug>

SparseMatrix::SparseMatrix(int row, int maxColumn):
  m_rowCount(row),
  m_maxColCount(maxColumn)
{
  m_rows = new int[m_rowCount];
  m_ncols = new int[m_rowCount];

  m_colIndexes = new int*[m_rowCount];
  m_values = new double*[m_rowCount];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_rowCount ; i++)
  {
    m_rows[i] = i;
    m_ncols[i] = 0;

    int *colIndexes = new int[m_maxColCount];

    for(int j = 0 ; j < m_maxColCount; j++)
    {
      colIndexes[j] = -1;
    }

    m_colIndexes[i] = colIndexes;
    m_values[i] = new double[m_maxColCount];
  }
}

SparseMatrix::~SparseMatrix()
{
  delete[] m_rows;
  delete[] m_ncols;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_rowCount ; i++)
  {
    delete[] m_colIndexes[i];
    delete[] m_values[i];
  }

  delete[] m_colIndexes;
  delete[] m_values;
}

const double &SparseMatrix::operator ()(int row, int column) const
{
  int *colIndexes = m_colIndexes[row];
  int ncols = m_ncols[row];

  for(int i = 0 ; i < ncols; i++)
  {
    if(colIndexes[i] == column)
    {
      return m_values[row][i];
    }
  }


  return m_zero;
}

double &SparseMatrix::operator ()(int row, int column)
{
  int *colIndexes = m_colIndexes[row];
  int ncols = m_ncols[row];

  for(int i = 0 ; i < ncols; i++)
  {
    if(colIndexes[i] == column)
    {
      return m_values[row][i];
    }
  }

  m_placeH = 0;
  return m_placeH;
}

void SparseMatrix::appendValue(int row, int column, double value)
{
  int ncols = m_ncols[row];
  m_values[row][ncols] = value;
  m_colIndexes[row][ncols] = column;
  m_ncols[row] = ncols + 1;
}

void SparseMatrix::setValue(int row, int column, double value)
{
  int *colIndexes = m_colIndexes[row];
  int ncols = m_ncols[row];

  for(int i = 0 ; i < ncols; i++)
  {
    if(colIndexes[i] == column)
    {
       m_values[row][i] = value;
       return;
    }
  }

  appendValue(row,column,value);
}

void SparseMatrix::getDataByRow(double data[]) const
{

  int c = 0;

  for(int i = 0 ; i < m_rowCount ; i++)
  {
    int ncols = m_ncols[i];

    for(int j = 0; j < ncols; j++)
    {
      data[c] = m_values[i][j];
      c++;
    }
  }
}

void SparseMatrix::getColumnIndexes(int indexes[]) const
{
  int c = 0;

  for(int i = 0 ; i < m_rowCount ; i++)
  {
    int ncols = m_ncols[i];

    for(int j = 0; j < ncols; j++)
    {
      indexes[c] = m_colIndexes[i][j];
      c++;
    }
  }
}

int SparseMatrix::getDataSize() const
{
  int count = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0; i < m_rowCount ; i++)
  {
#ifdef USE_OPENMP
#pragma omp atomic
#endif
    count += m_ncols[i];
  }

  return count;
}

int SparseMatrix::getColumnCount(int row) const
{
  return m_ncols[row];
}

Matrix<double> SparseMatrix::toMatrix() const
{
  Matrix<double> data(m_rowCount,m_rowCount);

  for(int i = 0; i < m_rowCount ; i++)
  {
    for(int j = 0; j < m_ncols[i]; j++)
    {
      data(i,m_colIndexes[i][j]) = m_values[i][j];
    }
  }

  return data;
}

int SparseMatrix::maxColumCount() const
{
  return m_maxColCount;
}

int SparseMatrix::numRows() const
{
  return m_rowCount;
}

void SparseMatrix::print() const
{
 toMatrix().printMatrix();
}



