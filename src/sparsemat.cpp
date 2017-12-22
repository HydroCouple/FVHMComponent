#include "stdafx.h"
#include "sparsemat.h"

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <QString>
#include <QDebug>

SparseMatrix::SparseMatrix(int ilower, int iupper, int columnCount, int maxColSize):
  m_rowCount(iupper-ilower + 1),
  m_columnCount(columnCount),
  m_maxColSize(maxColSize),
  m_ilower(ilower),
  m_iupper(iupper)
{
  m_ncols = new int[m_rowCount];
  m_colIndexes = new int*[m_rowCount];
  m_values = new double*[m_rowCount];

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = 0 ; i < m_rowCount ; i++)
  {
    m_ncols[i] = 0;

    int *colIndexes = new int[m_maxColSize];

    for(int j = 0 ; j < m_maxColSize; j++)
    {
      colIndexes[j] = -1;
    }

    m_colIndexes[i] = colIndexes;
    m_values[i] = new double[m_maxColSize];
  }
}

SparseMatrix::~SparseMatrix()
{
  delete[] m_ncols;

  for(int i = 0 ; i < m_rowCount ; i++)
  {
    delete[] m_colIndexes[i];
    delete[] m_values[i];
  }

  delete[] m_colIndexes;
  delete[] m_values;
}

int SparseMatrix::ilower() const
{
  return m_ilower;
}

int SparseMatrix::iupper() const
{
  return m_iupper;
}

int SparseMatrix::rowCount() const
{
  return m_rowCount;
}

int SparseMatrix::columnCount() const
{
  return m_columnCount;
}

int SparseMatrix::maxColumnSize() const
{
  return m_maxColSize;
}

const double &SparseMatrix::operator ()(int row, int column) const
{
  row = row - m_ilower;
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
  row = row - m_ilower;
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
  row = row - m_ilower;
  int ncols = m_ncols[row];
  m_values[row][ncols] = value;
  m_colIndexes[row][ncols] = column;
  m_ncols[row] = ncols + 1;
}

void SparseMatrix::setValue(int row, int column, double value)
{
  row = row - m_ilower;

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

int SparseMatrix::getDataSize(int ilower, int iupper) const
{
  int count = 0;

#ifdef USE_OPENMP
#pragma omp parallel for
#endif
  for(int i = ilower; i <= iupper ; i++)
  {
#ifdef USE_OPENMP
#pragma omp atomic
#endif
    count += m_ncols[i - m_ilower];
  }

  return count;
}

void SparseMatrix::getValuesByRows(double data[], int ilower, int iupper) const
{

  int c = 0;

  for(int i = ilower ; i <= iupper ; i++)
  {
    int rowIndex = i - m_ilower;
    int ncols = m_ncols[rowIndex];

    for(int j = 0; j < ncols; j++)
    {
      data[c] = m_values[rowIndex][j];
      c++;
    }
  }
}

void SparseMatrix::getColumnIndexes(int indexes[], int ilower, int iupper) const
{
  int c = 0;

  for(int i = ilower; i <= iupper ;i++)
  {
    int rowIndex = i - m_ilower;
    int ncols = m_ncols[rowIndex];

    for(int j = 0; j < ncols; j++)
    {
      indexes[c] = m_colIndexes[rowIndex][j];
      c++;
    }
  }
}

void SparseMatrix::getRowIndexes(int indexes[], int ilower, int iupper) const
{
  int c = 0;
  for(int i = ilower; i <= iupper; i++)
  {
    indexes[c] = i;
    c++;
  }
}

void SparseMatrix::getColsPerRow(int indexes[], int ilower, int iupper)  const
{
  int c = 0;
  for(int i = ilower; i <= iupper; i++)
  {
    indexes[c] = m_ncols[i - m_ilower];
    c++;
  }
}

int SparseMatrix::getColumnCount(int row) const
{
  return m_ncols[row - m_ilower];
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

void SparseMatrix::serializeRows(int ilower, int iupper, double values[], int &counter) const
{
  values[counter] = ilower; counter++;
  values[counter] = iupper; counter++;
  values[counter] = m_columnCount; counter++;
  values[counter] = m_maxColSize; counter++;

  for(int i = ilower; i <= iupper; i++)
  {
    int r = i - m_ilower;
    int ncols = m_ncols[r];
    values[counter] = ncols; counter++;

    for(int j = 0; j < ncols; j++)
    {
      values[counter] = m_colIndexes[r][j]; counter++;
      values[counter] = m_values[r][j]; counter++;
    }
  }
}

void SparseMatrix::deserialize(const double values[], SparseMatrix *&sparse, double *&b, double *&x, int &counter)
{

  int ilower = values[counter]; counter++;
  int iupper = values[counter]; counter++;
  int columnCount = values[counter]; counter++;
  int maxColumnSize = values[counter]; counter++;

  sparse = new SparseMatrix(ilower, iupper, columnCount, maxColumnSize);

  for(int i = ilower; i <= iupper ; i++)
  {
    int ncols = values[counter]; counter++;

    for(int j = 0; j < ncols; j++)
    {
      int colIndex = values[counter]; counter++;
      double value = values[counter]; counter++;
      sparse->appendValue(i,colIndex,value);
    }
  }

  b = new double[iupper - ilower + 1];

  for(int i = ilower; i <= iupper ; i++)
  {
    b[i-ilower] = values[counter]; counter++;
  }

  x = new double[iupper - ilower + 1];

  for(int i = ilower; i <= iupper ; i++)
  {
    x[i-ilower] = values[counter]; counter++;
  }
}

void SparseMatrix::print() const
{
 toMatrix().printMatrix();
}
