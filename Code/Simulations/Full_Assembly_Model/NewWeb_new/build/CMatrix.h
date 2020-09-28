// $Id: CMatrix.h 2500 2017-02-27 17:33:43Z axel $
//
// Coded by Alex Chirokov (there was no copyright notice).
// Enhanced by Axel Rossberg.
//
// Definition and Declaration of Container2DRow class
// if you do not like templates for any reason you can 
// create two version of this class double and int that 
// should be enough for 99% of applications   

#ifndef __CMATRIX__
#define __CMATRIX__

#include <string.h>

template <class T>
class Container2DRow;

template <class T, class ExtendedContainer2DRow = class Container2DRow<T> >
class CMatrix;

template <class T>
class Container2DRow
{
  friend class CMatrix<T, class ExtendedContainer2DRow >;
public:
  inline T& operator [] (int j);
  inline const T& operator [] (int j) const; 
  T **m_ppMatrix;
  int m_nXSize;
  mutable int i; //ROW (Y coord)
};
///Class container

template<class T> 
inline const T& Container2DRow<T>::operator [] (int j) const 
{
    ASSERT(j>=0 && j<m_nXSize); 
    return m_ppMatrix[i][j];
}

template<class T> 
inline T& Container2DRow<T>::operator [] (int j) 
{
    ASSERT(j>=0 && j<m_nXSize); 
    return m_ppMatrix[i][j];
}

/// Template for defining simple matrices.
template <class T, class ExtendedContainer2DRow >
class CMatrix  
{
public:
    //Helper class for [][] indexing, it is not neccesarily 
    // to agragated by CMatrix it could be just a friend
    ExtendedContainer2DRow row;

private:
    int m_nXSize;
    int m_nYSize;
    size_t m_nMemorySize;
    T **m_ppMatrix;

    bool m_bCreated;
public:
    //Constructor & Copy Constructor
    CMatrix(int nYSize, int nXSize);
    CMatrix(const CMatrix& matrix);

    //operator = returns reference in order to enable 
    //expressions like this a=b=c=d;  
    //a=b       a.operator=(b)
    //a=b+c     a.operator=(b.operator+(c));
    //a=b-c     a.operator=(b.operator-(c)); 
    CMatrix& operator= (const CMatrix& matrix);
    CMatrix& assign_filling(const CMatrix& matrix); //just as = but
						    //does not shrink
						    //size
    CMatrix  operator+ (const T& item);
    CMatrix  operator- (const T& item);
    CMatrix  operator* (const T& item);
    CMatrix  operator/ (const T& item);

    //Indexing //Y(row) X(col) 
    T& operator()(int i, int j) const;   // i - row
    //operator  [] returns object of type  Container2DRow
    //with have operator [] overloaded and know how to access 
    //matrix data 
#if 0//defined(DEBUGGING) && !defined(PARALLEL)
    inline ExtendedContainer2DRow & operator [] (int i);
    inline const    ExtendedContainer2DRow & operator [] (int i) const; 
#else
    inline T * operator [] (int i);
    inline const T * operator [] (int i) const; 
#endif

    //Helper functions, you can expand this section to do
    //LU decomposition, determinant evaluation and so on,  
    T SumAll() const;
    //Get Size
    int GetXSize() const;
    int GetYSize() const;
    T GetMinValue() const;
    T GetMaxValue() const;
    virtual ~CMatrix();
};
template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow>::CMatrix(int nYSize, int nXSize):row()
{
    m_bCreated = false;
    ASSERT(nXSize>=0 && nYSize>=0);


    m_nXSize = nXSize;
    m_nYSize = nYSize;
    m_nMemorySize = size_t(m_nYSize)*size_t(m_nXSize)*sizeof(T);

    m_ppMatrix    = new T* [m_nYSize ? m_nYSize : 1];
    m_ppMatrix[0] = new T  [m_nYSize*m_nXSize];

    for (int i=1; i<m_nYSize; i++)
        m_ppMatrix[i] = m_ppMatrix[0]+i*m_nXSize;

    memset(m_ppMatrix[0], 0, m_nMemorySize);
    m_bCreated = true;
    row.m_ppMatrix = m_ppMatrix;
    row.m_nXSize   = m_nXSize;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow>::CMatrix(const CMatrix& matrix):
  row(matrix.row)
{
    m_nXSize = matrix.m_nXSize;
    m_nYSize = matrix.m_nYSize;
    m_nMemorySize = m_nYSize*m_nXSize*sizeof(T);

    m_ppMatrix    = new T* [m_nYSize];
    ASSERT(m_ppMatrix!=NULL);

    m_ppMatrix[0] = new T  [m_nYSize*m_nXSize];
    ASSERT(m_ppMatrix[0]!=NULL);

    for (int i=1; i<m_nYSize; i++)
        m_ppMatrix[i] = m_ppMatrix[0]+i*m_nXSize;

    memcpy(m_ppMatrix[0],matrix.m_ppMatrix[0], m_nMemorySize);

    row.m_ppMatrix = m_ppMatrix;
    row.m_nXSize   = m_nXSize;

    m_bCreated = true;
}


template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow>& CMatrix<T,ExtendedContainer2DRow>::operator= (const CMatrix& matrix)
{
    if (this == &matrix) return *this;

    row=matrix.row;

    if(m_nXSize != matrix.m_nXSize || 
       m_nYSize != matrix.m_nYSize){
      // reallocate and organize memory:
      if (m_bCreated)
	{
	  delete [] m_ppMatrix[0];
	  delete [] m_ppMatrix;
	}
      m_nXSize = matrix.m_nXSize;
      m_nYSize = matrix.m_nYSize;
      m_nMemorySize = m_nYSize*m_nXSize*sizeof(T);

      m_ppMatrix    = new T* [m_nYSize];
      m_ppMatrix[0] = new T  [m_nYSize*m_nXSize];
      
      for (int i=1; i<m_nYSize; i++)
        m_ppMatrix[i] = m_ppMatrix[0]+i*m_nXSize;

      m_bCreated = true;
    }

    row.m_ppMatrix = m_ppMatrix;
    row.m_nXSize   = m_nXSize;

    memcpy(m_ppMatrix[0],matrix.m_ppMatrix[0], m_nMemorySize);
       
    return *this;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow>& CMatrix<T,ExtendedContainer2DRow>::assign_filling (const CMatrix& matrix)
{
    if (this == &matrix) return *this;

    row=matrix.row;

    if(m_nXSize < matrix.m_nXSize || 
       m_nYSize < matrix.m_nYSize){
      // reallocate and organize memory:
      if (m_bCreated)
	{
	  delete [] m_ppMatrix[0];
	  delete [] m_ppMatrix;
	}
      m_nXSize = matrix.m_nXSize;
      m_nYSize = matrix.m_nYSize;
      m_nMemorySize = m_nYSize*m_nXSize*sizeof(T);

      m_ppMatrix    = new T* [m_nYSize];
      m_ppMatrix[0] = new T  [m_nYSize*m_nXSize];
      
      for (int i=1; i<m_nYSize; i++)
        m_ppMatrix[i] = m_ppMatrix[0]+i*m_nXSize;

      m_bCreated = true;
    }

    row.m_ppMatrix = m_ppMatrix;
    row.m_nXSize   = m_nXSize;

    const size_t rowsize=size_t(matrix.m_nXSize)*sizeof(T);
    for (int i=0; i<matrix.m_nYSize; i++){
      memcpy(m_ppMatrix[i],matrix.m_ppMatrix[i],rowsize);
    }

    return *this;
}

template<class T, class ExtendedContainer2DRow >
T CMatrix<T,ExtendedContainer2DRow>::GetMinValue() const
{
    T minValue = m_ppMatrix[0][0];
    int i,j;

    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            if(m_ppMatrix[i][j]<minValue)
                minValue = m_ppMatrix[i][j];
        }
        return minValue;
}

template<class T, class ExtendedContainer2DRow >
T CMatrix<T,ExtendedContainer2DRow>::GetMaxValue() const
{
    T maxValue = m_ppMatrix[0][0];
    int i,j;

    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            if(m_ppMatrix[i][j]>maxValue)
                maxValue = m_ppMatrix[i][j];
        }
        return maxValue;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow> CMatrix<T,ExtendedContainer2DRow>::operator+ (const T& item)
{
    int i, j;

    CMatrix<T,ExtendedContainer2DRow> mtrx(m_nYSize, m_nXSize);
    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            mtrx.m_ppMatrix[i][j] = m_ppMatrix[i][j]+item ;
        }
        return mtrx;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow> CMatrix<T,ExtendedContainer2DRow>::operator- (const T& item)
{
    int i, j;

    CMatrix<T,ExtendedContainer2DRow> mtrx(m_nYSize, m_nXSize);
    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            mtrx.m_ppMatrix[i][j] = m_ppMatrix[i][j]-item ;
        }
        return mtrx;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow> CMatrix<T,ExtendedContainer2DRow>::operator* (const T& item)
{
    int i, j;

    CMatrix<T,ExtendedContainer2DRow> mtrx(m_nYSize, m_nXSize);
    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            mtrx.m_ppMatrix[i][j] = m_ppMatrix[i][j]*item ;
        }
        return mtrx;
}

template<class T, class ExtendedContainer2DRow >

CMatrix<T,ExtendedContainer2DRow> CMatrix<T,ExtendedContainer2DRow>::operator/ (const T& item)
{
    int i, j;

    CMatrix<T,ExtendedContainer2DRow> mtrx(m_nYSize, m_nXSize);
    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            mtrx.m_ppMatrix[i][j] = m_ppMatrix[i][j]*(1/item) ;
        }
        return mtrx;
}

template<class T, class ExtendedContainer2DRow >
CMatrix<T,ExtendedContainer2DRow>::~CMatrix()
{
    if (m_bCreated)
    {
        delete [] m_ppMatrix[0];
        delete [] m_ppMatrix;
    }
}

template<class T, class ExtendedContainer2DRow >
int CMatrix<T,ExtendedContainer2DRow>::GetXSize() const
{
    return m_nXSize;
}

template<class T, class ExtendedContainer2DRow >
T CMatrix<T,ExtendedContainer2DRow>::SumAll() const
{
    T sum = 0;
    int i, j;

    for (i=0; i<m_nYSize; i++)
        for (j=0; j<m_nXSize; j++)
        {
            sum += m_ppMatrix[i][j];
        }
        return sum;
}

template<class T, class ExtendedContainer2DRow >
int CMatrix<T,ExtendedContainer2DRow>::GetYSize() const
{
    return m_nYSize;
}
template<class T, class ExtendedContainer2DRow >        //Y(row) X(col)      
inline T& CMatrix<T,ExtendedContainer2DRow>::operator()(int i, int j) const
{
    ASSERT(i>=0 && i<m_nYSize &&
        j>=0 && j<m_nXSize);

    return m_ppMatrix[i][j];
}

#if 0//defined(DEBUGGING) && !defined(PARALLEL)
//Fancy Indexing
template<class T, class ExtendedContainer2DRow > 
inline ExtendedContainer2DRow & CMatrix<T,ExtendedContainer2DRow>::operator [] (int i)

{
    ASSERT(i>=0 && i<m_nYSize); 
    row.i = i;
    return row;
}

template<class T, class ExtendedContainer2DRow > 
inline const ExtendedContainer2DRow & CMatrix<T,ExtendedContainer2DRow>::operator [] (int i) const
{
    ASSERT(i>=0 && i<m_nYSize); 
    row.i = i;
    return row;
}
#else
//Fast Indexing, DOES NOT GO THROUGH ExtendedContainer2DRow
template<class T, class ExtendedContainer2DRow > 
inline T *  CMatrix<T,ExtendedContainer2DRow>::operator [] (int i)

{
    return m_ppMatrix[i];
}

template<class T, class ExtendedContainer2DRow > 
inline const T *  CMatrix<T,ExtendedContainer2DRow>::operator [] (int i) const
{
    return m_ppMatrix[i];
}
#endif // DEBUGGING

#endif // __CMATRIX__
