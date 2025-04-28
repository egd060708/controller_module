// $Id: Array.hh 249 2008-11-20 09:58:23Z schaerf $
// This file is part of EasyLocalpp: a C++ Object-Oriented framework
// aimed at easing the development of Local Search algorithms.
// Copyright (C) 2001--2008 Andrea Schaerf, Luca Di Gaspero. 
//
// This software may be modified and distributed under the terms
// of the MIT license.  See the LICENSE file for details.

#if !defined(_ARRAY_HH)
#define _ARRAY_HH

#include <set>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

namespace quadprogpp {

enum MType { DIAG };

template <typename T>
class qpVector
{
public: 
  qpVector(); 
  qpVector(const unsigned int n);  
  qpVector(const T& a, const unsigned int n); //initialize to constant value 
  qpVector(const T* a, const unsigned int n); // Initialize to array 
  qpVector(const qpVector &rhs); // copy constructor 
  ~qpVector(); // destructor
	
  inline void set(const T* a, const unsigned int n);
  qpVector<T> extract(const std::set<unsigned int>& indexes) const;
  inline T& operator[](const unsigned int& i); //i-th element 
  inline const T& operator[](const unsigned int& i) const; 
	
  inline unsigned int size() const;
  inline void resize(const unsigned int n);
  inline void resize(const T& a, const unsigned int n);
	
  qpVector<T>& operator=(const qpVector<T>& rhs); //assignment 
  qpVector<T>& operator=(const T& a); //assign a to every element 
  inline qpVector<T>& operator+=(const qpVector<T>& rhs);
  inline qpVector<T>& operator-=(const qpVector<T>& rhs);
  inline qpVector<T>& operator*=(const qpVector<T>& rhs);
  inline qpVector<T>& operator/=(const qpVector<T>& rhs);
  inline qpVector<T>& operator^=(const qpVector<T>& rhs);
  inline qpVector<T>& operator+=(const T& a);
  inline qpVector<T>& operator-=(const T& a);
  inline qpVector<T>& operator*=(const T& a);
  inline qpVector<T>& operator/=(const T& a);
  inline qpVector<T>& operator^=(const T& a);
private: 
  unsigned int n; // size of array. upper index is n-1 
  T* v; // storage for data
}; 

template <typename T> 
qpVector<T>::qpVector() 
  : n(0), v(0) 
{} 

template <typename T> 
qpVector<T>::qpVector(const unsigned int n) 
  : v(new T[n]) 
{
  this->n = n;
} 

template <typename T> 
qpVector<T>::qpVector(const T& a, const unsigned int n) 
  : v(new T[n])
{ 
  this->n = n;
  for (unsigned int i = 0; i < n; i++) 
    v[i] = a; 
} 

template <typename T> 
qpVector<T>::qpVector(const T* a, const unsigned int n) 
  : v(new T[n])
{ 
  this->n = n;
  for (unsigned int i = 0; i < n; i++) 
    v[i] = *a++; 
} 

template <typename T> 
qpVector<T>::qpVector(const qpVector<T>& rhs) 
  : v(new T[rhs.n])
{ 
  this->n = rhs.n;
  for (unsigned int	i = 0; i < n; i++) 
    v[i] = rhs[i]; 
} 

template <typename T> 
qpVector<T>::~qpVector() 
{ 
  if (v != 0) 
    delete[] (v); 
} 

template <typename T> 
void qpVector<T>::resize(const unsigned int n) 
{
  if (n == this->n)
    return;
  if (v != 0) 
    delete[] (v); 
  v = new T[n];
  this->n = n;
} 

template <typename T> 
void qpVector<T>::resize(const T& a, const unsigned int n) 
{
  resize(n);
  for (unsigned int i = 0; i < n; i++)
    v[i] = a;
} 


template <typename T> 
inline qpVector<T>& qpVector<T>::operator=(const qpVector<T>& rhs) 
// postcondition: normal assignment via copying has been performed; 
// if vector and rhs were different sizes, vector 
// has been resized to match the size of rhs 
{ 
  if (this != &rhs) 
    { 
      resize(rhs.n);
      for (unsigned int i = 0; i < n; i++) 
	v[i] = rhs[i]; 
    } 
  return *this; 
} 

template <typename T> 
inline qpVector<T> & qpVector<T>::operator=(const T& a) //assign a to every element 
{ 
  for (unsigned int i = 0; i < n; i++) 
    v[i] = a; 
  return *this; 
} 

template <typename T> 
inline T & qpVector<T>::operator[](const unsigned int& i) //subscripting 
{ 
  return v[i]; 
}

template <typename T>
inline const T& qpVector<T>::operator[](const unsigned int& i) const //subscripting 
{ 
  return v[i]; 
} 

template <typename T> 
inline unsigned int qpVector<T>::size() const 
{ 
  return n; 
}

template <typename T> 
inline void qpVector<T>::set(const T* a, unsigned int n) 
{ 
  resize(n);
  for (unsigned int i = 0; i < n; i++) 
    v[i] = a[i]; 
} 

template <typename T> 
inline qpVector<T> qpVector<T>::extract(const std::set<unsigned int>& indexes) const
{
  qpVector<T> tmp(indexes.size());
  unsigned int i = 0;
	
  for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
    {
      if (*el >= n)
	throw std::logic_error("Error extracting subvector: the indexes are out of vector bounds");
      tmp[i++] = v[*el];
    }
	
  return tmp;
}

template <typename T> 
inline qpVector<T>& qpVector<T>::operator+=(const qpVector<T>& rhs)
{
  if (this->size() != rhs.size())
    throw std::logic_error("Operator+=: vectors have different sizes");
  for (unsigned int i = 0; i < n; i++)
    v[i] += rhs[i];
	
  return *this;
}


template <typename T> 
inline qpVector<T>& qpVector<T>::operator+=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    v[i] += a;
	
  return *this;
}

template <typename T>
inline qpVector<T> operator+(const qpVector<T>& rhs)
{
  return rhs;
}

template <typename T>
inline qpVector<T> operator+(const qpVector<T>& lhs, const qpVector<T>& rhs)
{
  if (lhs.size() != rhs.size())
    throw std::logic_error("Operator+: vectors have different sizes");
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] + rhs[i];
	
  return tmp;
}

template <typename T>
inline qpVector<T> operator+(const qpVector<T>& lhs, const T& a)
{
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] + a;
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator+(const T& a, const qpVector<T>& rhs)
{
  qpVector<T> tmp(rhs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    tmp[i] = a + rhs[i];
		
  return tmp;
}

template <typename T> 
inline qpVector<T>& qpVector<T>::operator-=(const qpVector<T>& rhs)
{
  if (this->size() != rhs.size())
    throw std::logic_error("Operator-=: vectors have different sizes");
  for (unsigned int i = 0; i < n; i++)
    v[i] -= rhs[i];
	
  return *this;
}


template <typename T> 
inline qpVector<T>& qpVector<T>::operator-=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    v[i] -= a;
	
  return *this;
}

template <typename T>
inline qpVector<T> operator-(const qpVector<T>& rhs)
{
  return (T)(-1) * rhs;
}

template <typename T>
inline qpVector<T> operator-(const qpVector<T>& lhs, const qpVector<T>& rhs)
{
  if (lhs.size() != rhs.size())
    throw std::logic_error("Operator-: vectors have different sizes");
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] - rhs[i];
	
  return tmp;
}

template <typename T>
inline qpVector<T> operator-(const qpVector<T>& lhs, const T& a)
{
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] - a;
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator-(const T& a, const qpVector<T>& rhs)
{
  qpVector<T> tmp(rhs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    tmp[i] = a - rhs[i];
		
  return tmp;
}

template <typename T> 
inline qpVector<T>& qpVector<T>::operator*=(const qpVector<T>& rhs)
{
  if (this->size() != rhs.size())
    throw std::logic_error("Operator*=: vectors have different sizes");
  for (unsigned int i = 0; i < n; i++)
    v[i] *= rhs[i];
	
  return *this;
}


template <typename T> 
inline qpVector<T>& qpVector<T>::operator*=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    v[i] *= a;
	
  return *this;
}

template <typename T>
inline qpVector<T> operator*(const qpVector<T>& lhs, const qpVector<T>& rhs)
{
  if (lhs.size() != rhs.size())
    throw std::logic_error("Operator*: vectors have different sizes");
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] * rhs[i];
	
  return tmp;
}

template <typename T>
inline qpVector<T> operator*(const qpVector<T>& lhs, const T& a)
{
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] * a;
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator*(const T& a, const qpVector<T>& rhs)
{
  qpVector<T> tmp(rhs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    tmp[i] = a * rhs[i];
		
  return tmp;
}

template <typename T> 
inline qpVector<T>& qpVector<T>::operator/=(const qpVector<T>& rhs)
{
  if (this->size() != rhs.size())
    throw std::logic_error("Operator/=: vectors have different sizes");
  for (unsigned int i = 0; i < n; i++)
    v[i] /= rhs[i];
	
  return *this;
}


template <typename T> 
inline qpVector<T>& qpVector<T>::operator/=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    v[i] /= a;
	
  return *this;
}

template <typename T>
inline qpVector<T> operator/(const qpVector<T>& lhs, const qpVector<T>& rhs)
{
  if (lhs.size() != rhs.size())
    throw std::logic_error("Operator/: vectors have different sizes");
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] / rhs[i];
	
  return tmp;
}

template <typename T>
inline qpVector<T> operator/(const qpVector<T>& lhs, const T& a)
{
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = lhs[i] / a;
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator/(const T& a, const qpVector<T>& rhs)
{
  qpVector<T> tmp(rhs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    tmp[i] = a / rhs[i];
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator^(const qpVector<T>& lhs, const qpVector<T>& rhs)
{
  if (lhs.size() != rhs.size())
    throw std::logic_error("Operator^: vectors have different sizes");
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = pow(lhs[i], rhs[i]);
	
  return tmp;
}

template <typename T>
inline qpVector<T> operator^(const qpVector<T>& lhs, const T& a)
{
  qpVector<T> tmp(lhs.size());
  for (unsigned int i = 0; i < lhs.size(); i++)
    tmp[i] = pow(lhs[i], a);
		
  return tmp;
}

template <typename T>
inline qpVector<T> operator^(const T& a, const qpVector<T>& rhs)
{
  qpVector<T> tmp(rhs.size());
  for (unsigned int i = 0; i < rhs.size(); i++)
    tmp[i] = pow(a, rhs[i]);
		
  return tmp;
}

template <typename T>
inline qpVector<T>& qpVector<T>::operator^=(const qpVector<T>& rhs)
{
  if (this->size() != rhs.size())
    throw std::logic_error("Operator^=: vectors have different sizes");
  for (unsigned int i = 0; i < n; i++)
    v[i] = pow(v[i], rhs[i]);
		
  return *this;
}

template <typename T>
inline qpVector<T>& qpVector<T>::operator^=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    v[i] = pow(v[i], a);
		
  return *this;
}

template <typename T>
inline bool operator==(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] != w[i])
      return false;
  return true;
}

template <typename T>
inline bool operator!=(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] != w[i])
      return true;
  return false;
}

template <typename T>
inline bool operator<(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] >= w[i])
      return false;
  return true;
}

template <typename T>
inline bool operator<=(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] > w[i])
      return false;
  return true;
}

template <typename T>
inline bool operator>(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] <= w[i])
      return false;
  return true;
}

template <typename T>
inline bool operator>=(const qpVector<T>& v, const qpVector<T>& w)
{
  if (v.size() != w.size())
    throw std::logic_error("Vectors of different size are not confrontable");
  for (unsigned i = 0; i < v.size(); i++)
    if (v[i] < w[i])
      return false;
  return true;
}

/**
   Input/Output 
*/
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const qpVector<T>& v)
{
  os << std::endl << v.size() << std::endl;
  for (unsigned int i = 0; i < v.size() - 1; i++)
    os << std::setw(20) << std::setprecision(16) << v[i] << ", ";
  os << std::setw(20) << std::setprecision(16) << v[v.size() - 1] << std::endl;
	
  return os;
}

template <typename T>
std::istream& operator>>(std::istream& is, qpVector<T>& v)
{
  int elements;
  char comma;
  is >> elements;
  v.resize(elements);
  for (unsigned int i = 0; i < elements; i++)
    is >> v[i] >> comma;
	
  return is;
}

/**
   Index utilities
*/

std::set<unsigned int> seq(unsigned int s, unsigned int e);

std::set<unsigned int> singleton(unsigned int i);

template <typename T>
class CanonicalBaseVector : public qpVector<T>
{
public:
  CanonicalBaseVector(unsigned int i, unsigned int n);
  inline void reset(unsigned int i);
private:
  unsigned int e;
};

template <typename T>
CanonicalBaseVector<T>::CanonicalBaseVector(unsigned int i, unsigned int n)
  : qpVector<T>((T)0, n), e(i)
{ (*this)[e] = (T)1; }

template <typename T>
inline void CanonicalBaseVector<T>::reset(unsigned int i)
{ 
  (*this)[e] = (T)0; 
  e = i; 
  (*this)[e] = (T)1;
}

#include <stdexcept>

template <typename T>
inline T sum(const qpVector<T>& v)
{
  T tmp = (T)0;
  for (unsigned int i = 0; i < v.size(); i++)
    tmp += v[i];
	
  return tmp;
}

template <typename T>
inline T prod(const qpVector<T>& v)
{
  T tmp = (T)1;
  for (unsigned int i = 0; i < v.size(); i++)
    tmp *= v[i];
	
  return tmp;
}

template <typename T>
inline T mean(const qpVector<T>& v)
{
  T sum = (T)0;
  for (unsigned int i = 0; i < v.size(); i++)
    sum += v[i];
  return sum / v.size();
}

template <typename T>
inline T median(const qpVector<T>& v)
{
  qpVector<T> tmp = sort(v);
  if (v.size() % 2 == 1) // it is an odd-sized vector
    return tmp[v.size() / 2];
  else
    return 0.5 * (tmp[v.size() / 2 - 1] + tmp[v.size() / 2]);
}

template <typename T>
inline T stdev(const qpVector<T>& v, bool sample_correction = false)
{
  return sqrt(var(v, sample_correction));
}

template <typename T>
inline T var(const qpVector<T>& v, bool sample_correction = false)
{
  T sum = (T)0, ssum = (T)0;
  unsigned int n = v.size();
  for (unsigned int i = 0; i < n; i++)
    {	
      sum += v[i];
      ssum += (v[i] * v[i]);
    }
  if (!sample_correction)
    return (ssum / n) - (sum / n) * (sum / n);
  else
    return n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
}

//template <typename T>
//inline T max(const qpVector<T>& v)
//{
//  T value = v[0];
//  for (unsigned int i = 1; i < v.size(); i++)
//    value = (std::max)(v[i], value);
//	
//  return value;
//}

//template <typename T>
//inline T min(const qpVector<T>& v)
//{
//  T value = v[0];
//  for (unsigned int i = 1; i < v.size(); i++)
//    value = (std::min)(v[i], value);
//	
//  return value;
//}

template <typename T>
inline unsigned int index_max(const qpVector<T>& v)
{
  unsigned int max = 0;
  for (unsigned int i = 1; i < v.size(); i++)
    if (v[i] > v[max])
      max = i;
	
  return max;
}

template <typename T>
inline unsigned int index_min(const qpVector<T>& v)
{
  unsigned int min = 0;
  for (unsigned int i = 1; i < v.size(); i++)
    if (v[i] < v[min])
      min = i;
	
  return min;
}


template <typename T>
inline T dot_prod(const qpVector<T>& a, const qpVector<T>& b)
{
  T sum = (T)0;
  if (a.size() != b.size())
    throw std::logic_error("Dotprod error: the vectors are not the same size");
  for (unsigned int i = 0; i < a.size(); i++)
    sum += a[i] * b[i];
	
  return sum;
}

/**
   Single element mathematical functions
*/

template <typename T>
inline qpVector<T> exp(const qpVector<T>& v)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = exp(v[i]);
	
  return tmp;
}

template <typename T>
inline qpVector<T> log(const qpVector<T>& v)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = log(v[i]);
	
  return tmp;
}

template <typename T>
inline qpVector<T> vec_sqrt(const qpVector<T>& v)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = sqrt(v[i]);
	
  return tmp;
}

template <typename T>
inline qpVector<T> pow(const qpVector<T>& v, double a)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = pow(v[i], a);
	
  return tmp;
}

template <typename T>
inline qpVector<T> abs(const qpVector<T>& v)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = (T)fabs(v[i]);
	
  return tmp;
}

template <typename T>
inline qpVector<T> sign(const qpVector<T>& v)
{
  qpVector<T> tmp(v.size());
  for (unsigned int i = 0; i < v.size(); i++)
    tmp[i] = v[i] > 0 ? +1 : v[i] == 0 ? 0 : -1;
	
  return tmp;
}

template <typename T>
inline unsigned int partition(qpVector<T>& v, unsigned int begin, unsigned int end)
{
  unsigned int i = begin + 1, j = begin + 1;
  T pivot = v[begin];
  while (j <= end) 
    {
      if (v[j] < pivot) {
	std::swap(v[i], v[j]);
	i++;
      }
      j++;
    }
  v[begin] = v[i - 1];
  v[i - 1] = pivot;
  return i - 2;
}
	

template <typename T>
inline void quicksort(qpVector<T>& v, unsigned int begin, unsigned int end)
{
  if (end > begin)
    {
      unsigned int index = partition(v, begin, end);
      quicksort(v, begin, index);
      quicksort(v, index + 2, end);
    }
}

template <typename T>
inline qpVector<T> sort(const qpVector<T>& v)
{
  qpVector<T> tmp(v);
  
  quicksort<T>(tmp, 0, tmp.size() - 1);
  
  return tmp;
}

template <typename T>
inline qpVector<double> rank(const qpVector<T>& v)
{
  qpVector<T> tmp(v);
  qpVector<double> tmp_rank(0.0, v.size());	
	
  for (unsigned int i = 0; i < tmp.size(); i++)
    {
      unsigned int smaller = 0, equal = 0;
      for (unsigned int j = 0; j < tmp.size(); j++)
	if (i == j)
	  continue;
	else
	  if (tmp[j] < tmp[i])
	    smaller++;
	  else if (tmp[j] == tmp[i])
	    equal++;
      tmp_rank[i] = smaller + 1;
      if (equal > 0)
	{
	  for (unsigned int j = 1; j <= equal; j++)
	    tmp_rank[i] += smaller + 1 + j;
	  tmp_rank[i] /= (double)(equal + 1);
	}
    }
	
  return tmp_rank;
}

//enum MType { DIAG };

template <typename T>
class qpMatrix 
{
public:
  qpMatrix(); // Default constructor
  qpMatrix(const unsigned int n, const unsigned int m); // Construct a n x m matrix
  qpMatrix(const T& a, const unsigned int n, const unsigned int m); // Initialize the content to constant a
  qpMatrix(MType t, const T& a, const T& o, const unsigned int n, const unsigned int m);
  qpMatrix(MType t, const qpVector<T>& v, const T& o, const unsigned int n, const unsigned int m);
  qpMatrix(const T* a, const unsigned int n, const unsigned int m); // Initialize to array 
  qpMatrix(const qpMatrix<T>& rhs); // Copy constructor
  ~qpMatrix(); // destructor
	
  inline T* operator[](const unsigned int& i) { return v[i]; } // Subscripting: row i
  inline const T* operator[](const unsigned int& i) const { return v[i]; }; // const subsctipting
	
  inline void resize(const unsigned int n, const unsigned int m);
  inline void resize(const T& a, const unsigned int n, const unsigned int m);
	
	
  inline qpVector<T> extractRow(const unsigned int i) const; 
  inline qpVector<T> extractColumn(const unsigned int j) const;
  inline qpVector<T> extractDiag() const;
  inline qpMatrix<T> extractRows(const std::set<unsigned int>& indexes) const;
  inline qpMatrix<T> extractColumns(const std::set<unsigned int>& indexes) const;
  inline qpMatrix<T> extract(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes) const;
	
  inline void set(const T* a, unsigned int n, unsigned int m);
  inline void set(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes, const qpMatrix<T>& m);
  inline void setRow(const unsigned int index, const qpVector<T>& v);
  inline void setRow(const unsigned int index, const qpMatrix<T>& v);
  inline void setRows(const std::set<unsigned int>& indexes, const qpMatrix<T>& m);
  inline void setColumn(const unsigned int index, const qpVector<T>& v);
  inline void setColumn(const unsigned int index, const qpMatrix<T>& v);
  inline void setColumns(const std::set<unsigned int>& indexes, const qpMatrix<T>& m);
	
	
  inline unsigned int nrows() const { return n; } // number of rows
  inline unsigned int ncols() const { return m; } // number of columns
	
  inline qpMatrix<T>& operator=(const qpMatrix<T>& rhs); // Assignment operator
  inline qpMatrix<T>& operator=(const T& a); // Assign to every element value a
  inline qpMatrix<T>& operator+=(const qpMatrix<T>& rhs);
  inline qpMatrix<T>& operator-=(const qpMatrix<T>& rhs);
  inline qpMatrix<T>& operator*=(const qpMatrix<T>& rhs);
  inline qpMatrix<T>& operator/=(const qpMatrix<T>& rhs);
  inline qpMatrix<T>& operator^=(const qpMatrix<T>& rhs);
  inline qpMatrix<T>& operator+=(const T& a);
  inline qpMatrix<T>& operator-=(const T& a);
  inline qpMatrix<T>& operator*=(const T& a);
  inline qpMatrix<T>& operator/=(const T& a);
  inline qpMatrix<T>& operator^=(const T& a);
  inline operator qpVector<T>();
private:
  unsigned int n; // number of rows
  unsigned int m; // number of columns
  T **v; // storage for data
};

template <typename T>
qpMatrix<T>::qpMatrix() 
  : n(0), m(0), v(0)
{}

template <typename T>
qpMatrix<T>::qpMatrix(unsigned int n, unsigned int m)
  : v(new T*[n])
{
  this->n = n; this->m = m;
  v[0] = new T[m * n];
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
}

template <typename T>
qpMatrix<T>::qpMatrix(const T& a, unsigned int n, unsigned int m)
  : v(new T*[n])
{
  this->n = n; this->m = m;
  v[0] = new T[m * n];
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = a;
}

template <class T> 
qpMatrix<T>::qpMatrix(const T* a, unsigned int n, unsigned int m) 
  : v(new T*[n])
{ 
  this->n = n; this->m = m;
  v[0] = new T[m * n]; 
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m; 
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = *a++; 
} 

template <class T> 
qpMatrix<T>::qpMatrix(MType t, const T& a, const T& o, unsigned int n, unsigned int m) 
  : v(new T*[n])
{ 
  this->n = n; this->m = m;
  v[0] = new T[m * n]; 
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m; 
  switch (t)
    {
    case DIAG:
      for (unsigned int i = 0; i < n; i++)
	for (unsigned int j = 0; j < m; j++)
	  if (i != j)
	    v[i][j] = o; 
	  else
	    v[i][j] = a;
      break;
    default:
      throw std::logic_error("qpMatrix type not supported");
    }
} 

template <class T> 
qpMatrix<T>::qpMatrix(MType t, const qpVector<T>& a, const T& o, unsigned int n, unsigned int m) 
  : v(new T*[n])
{ 
  this->n = n; this->m = m;
  v[0] = new T[m * n]; 
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m; 
  switch (t)
    {
    case DIAG:
      for (unsigned int i = 0; i < n; i++)
	for (unsigned int j = 0; j < m; j++)
	  if (i != j)
	    v[i][j] = o; 
	  else
	    v[i][j] = a[i];
      break;
    default:
      throw std::logic_error("qpMatrix type not supported");
    }
} 

template <typename T>
qpMatrix<T>::qpMatrix(const qpMatrix<T>& rhs)
  : v(new T*[rhs.n])
{
  n = rhs.n; m = rhs.m;
  v[0] = new T[m * n]; 
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = rhs[i][j];
}

template <typename T> 
qpMatrix<T>::~qpMatrix() 
{ 
  if (v != 0) { 
    delete[] (v[0]); 
    delete[] (v); 
  } 
}
				
template <typename T> 
inline qpMatrix<T>& qpMatrix<T>::operator=(const qpMatrix<T> &rhs) 
// postcondition: normal assignment via copying has been performed; 
// if matrix and rhs were different sizes, matrix 
// has been resized to match the size of rhs 
{ 
  if (this != &rhs) 
    {
      resize(rhs.n, rhs.m);
      for (unsigned int i = 0; i < n; i++)
	for (unsigned int j = 0; j < m; j++)
	  v[i][j] = rhs[i][j]; 
    } 
  return *this; 
} 

template <typename T> 
inline qpMatrix<T>& qpMatrix<T>::operator=(const T& a) // assign a to every element 
{ 
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = a; 
  return *this; 
} 


template <typename T> 
inline void qpMatrix<T>::resize(const unsigned int n, const unsigned int m) 
{
  if (n == this->n && m == this->m)
    return;
  if (v != 0) 
    { 
      delete[] (v[0]); 
      delete[] (v); 
    } 
  this->n = n; this->m = m;
  v = new T*[n]; 
  v[0] = new T[m * n];  
  for (unsigned int i = 1; i < n; i++)
    v[i] = v[i - 1] + m;
} 

template <typename T> 
inline void qpMatrix<T>::resize(const T& a, const unsigned int n, const unsigned int m) 
{
  resize(n, m);
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = a;
} 



template <typename T> 
inline qpVector<T> qpMatrix<T>::extractRow(const unsigned int i) const
{
  if (i >= n)
    throw std::logic_error("Error in extractRow: trying to extract a row out of matrix bounds");
  qpVector<T> tmp(v[i], m);
	
  return tmp;
}

template <typename T> 
inline qpVector<T> qpMatrix<T>::extractColumn(const unsigned int j) const
{
  if (j >= m)
    throw std::logic_error("Error in extractRow: trying to extract a row out of matrix bounds");
  qpVector<T> tmp(n);
	
  for (unsigned int i = 0; i < n; i++)
    tmp[i] = v[i][j];
	
  return tmp;
}

template <typename T>
inline qpVector<T> qpMatrix<T>::extractDiag() const
{
  unsigned int d = ((std::min))(n, m);
  
  qpVector<T> tmp(d);
	
  for (unsigned int i = 0; i < d; i++)
    tmp[i] = v[i][i];
	
  return tmp;
	
}

template <typename T> 
inline qpMatrix<T> qpMatrix<T>::extractRows(const std::set<unsigned int>& indexes) const
{
  qpMatrix<T> tmp(indexes.size(), m);
  unsigned int i = 0;
	
  for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
    {
      for (unsigned int j = 0; j < m; j++)
	{
	  if (*el >= n)
	    throw std::logic_error("Error extracting rows: the indexes are out of matrix bounds");
	  tmp[i][j] = v[*el][j];
	}
      i++;
    }
	
  return tmp;
}

template <typename T> 
inline qpMatrix<T> qpMatrix<T>::extractColumns(const std::set<unsigned int>& indexes) const
{
  qpMatrix<T> tmp(n, indexes.size());
  unsigned int j = 0;
	
  for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
    {
      for (unsigned int i = 0; i < n; i++)
	{
	  if (*el >= m)
	    throw std::logic_error("Error extracting columns: the indexes are out of matrix bounds");
	  tmp[i][j] = v[i][*el];
	}
      j++;
    }
	
  return tmp;
}

template <typename T> 
inline qpMatrix<T> qpMatrix<T>::extract(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes) const
{
  qpMatrix<T> tmp(r_indexes.size(), c_indexes.size());
  unsigned int i = 0, j;
	
  for (std::set<unsigned int>::const_iterator r_el = r_indexes.begin(); r_el != r_indexes.end(); r_el++)
    {
      if (*r_el >= n)
	throw std::logic_error("Error extracting submatrix: the indexes are out of matrix bounds");
      j = 0;
      for (std::set<unsigned int>::const_iterator c_el = c_indexes.begin(); c_el != c_indexes.end(); c_el++)
	{
	  if (*c_el >= m)
	    throw std::logic_error("Error extracting rows: the indexes are out of matrix bounds");
	  tmp[i][j] = v[*r_el][*c_el];
	  j++;
	}
      i++;
    }
	
  return tmp;
}

template <typename T> 
inline void qpMatrix<T>::setRow(unsigned int i, const qpVector<T>& a)
{	
  if (i >= n)
    throw std::logic_error("Error in setRow: trying to set a row out of matrix bounds");
  if (this->m != a.size())
    throw std::logic_error("Error setting matrix row: ranges are not compatible");
  for (unsigned int j = 0; j < ncols(); j++)
    v[i][j] = a[j];
}

template <typename T> 
inline void qpMatrix<T>::setRow(unsigned int i, const qpMatrix<T>& a)
{	
  if (i >= n)
    throw std::logic_error("Error in setRow: trying to set a row out of matrix bounds");
  if (this->m != a.ncols())
    throw std::logic_error("Error setting matrix column: ranges are not compatible");
  if (a.nrows() != 1)
    throw std::logic_error("Error setting matrix column with a non-row matrix");
  for (unsigned int j = 0; j < ncols(); j++)
    v[i][j] = a[0][j];
}

template <typename T> 
inline void qpMatrix<T>::setRows(const std::set<unsigned int>& indexes, const qpMatrix<T>& m)
{
  unsigned int i = 0;
	
  if (indexes.size() != m.nrows() || this->m != m.ncols())
    throw std::logic_error("Error setting matrix rows: ranges are not compatible");
  for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
    {
      for (unsigned int j = 0; j < ncols(); j++)
	{
	  if (*el >= n)
	    throw std::logic_error("Error in setRows: trying to set a row out of matrix bounds");
	  v[*el][j] = m[i][j];
	}
      i++;
    }
}

template <typename T> 
inline void qpMatrix<T>::setColumn(unsigned int j, const qpVector<T>& a)
{	
  if (j >= m)
    throw std::logic_error("Error in setColumn: trying to set a column out of matrix bounds");
  if (this->n != a.size())
    throw std::logic_error("Error setting matrix column: ranges are not compatible");
  for (unsigned int i = 0; i < nrows(); i++)
    v[i][j] = a[i];
}

template <typename T> 
inline void qpMatrix<T>::setColumn(unsigned int j, const qpMatrix<T>& a)
{	
  if (j >= m)
    throw std::logic_error("Error in setColumn: trying to set a column out of matrix bounds");
  if (this->n != a.nrows())
    throw std::logic_error("Error setting matrix column: ranges are not compatible");
  if (a.ncols() != 1)
    throw std::logic_error("Error setting matrix column with a non-column matrix");
  for (unsigned int i = 0; i < nrows(); i++)
    v[i][j] = a[i][0];
}


template <typename T> 
inline void qpMatrix<T>::setColumns(const std::set<unsigned int>& indexes, const qpMatrix<T>& a)
{
  unsigned int j = 0;
	
  if (indexes.size() != a.ncols() || this->n != a.nrows())
    throw std::logic_error("Error setting matrix columns: ranges are not compatible");
  for (std::set<unsigned int>::const_iterator el = indexes.begin(); el != indexes.end(); el++)
    {
      for (unsigned int i = 0; i < nrows(); i++)
	{
	  if (*el >= m)
	    throw std::logic_error("Error in setColumns: trying to set a column out of matrix bounds");
	  v[i][*el] = a[i][j];
	}
      j++;
    }
}

template <typename T> 
inline void qpMatrix<T>::set(const std::set<unsigned int>& r_indexes, const std::set<unsigned int>& c_indexes, const qpMatrix<T>& a)
{
  unsigned int i = 0, j;
  if (c_indexes.size() != a.ncols() || r_indexes.size() != a.nrows())
    throw std::logic_error("Error setting matrix elements: ranges are not compatible");
	
  for (std::set<unsigned int>::const_iterator r_el = r_indexes.begin(); r_el != r_indexes.end(); r_el++)
    {
      if (*r_el >= n)
	throw std::logic_error("Error in set: trying to set a row out of matrix bounds");
      j = 0;
      for (std::set<unsigned int>::const_iterator c_el = c_indexes.begin(); c_el != c_indexes.end(); c_el++)
	{
	  if (*c_el >= m)
	    throw std::logic_error("Error in set: trying to set a column out of matrix bounds");
	  v[*r_el][*c_el] = a[i][j];
	  j++;
	}
      i++;
    }
}

template <typename T> 
inline void qpMatrix<T>::set(const T* a, unsigned int n, unsigned int m)
{
  if (this->n != n || this->m != m)
    resize(n, m);
  unsigned int k = 0;
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = a[k++];
}


template <typename T>
qpMatrix<T> operator+(const qpMatrix<T>& rhs)
{
  return rhs;
}

template <typename T>
qpMatrix<T> operator+(const qpMatrix<T>& lhs, const qpMatrix<T>& rhs)
{
  if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
    throw std::logic_error("Operator+: matrices have different sizes");
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] + rhs[i][j];
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator+(const qpMatrix<T>& lhs, const T& a)
{
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] + a;
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator+(const T& a, const qpMatrix<T>& rhs)
{
  qpMatrix<T> tmp(rhs.nrows(), rhs.ncols());
  for (unsigned int i = 0; i < rhs.nrows(); i++)
    for (unsigned int j = 0; j < rhs.ncols(); j++)
      tmp[i][j] = a + rhs[i][j];
	
  return tmp;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator+=(const qpMatrix<T>& rhs)
{
  if (m != rhs.ncols() || n != rhs.nrows())
    throw std::logic_error("Operator+=: matrices have different sizes");
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] += rhs[i][j];
	
  return *this;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator+=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] += a;
	
  return *this;
}

template <typename T>
qpMatrix<T> operator-(const qpMatrix<T>& rhs)
{	
  return (T)(-1) * rhs;
}

template <typename T>
qpMatrix<T> operator-(const qpMatrix<T>& lhs, const qpMatrix<T>& rhs)
{
  if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
    throw std::logic_error("Operator-: matrices have different sizes");
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] - rhs[i][j];
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator-(const qpMatrix<T>& lhs, const T& a)
{
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] - a;
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator-(const T& a, const qpMatrix<T>& rhs)
{
  qpMatrix<T> tmp(rhs.nrows(), rhs.ncols());
  for (unsigned int i = 0; i < rhs.nrows(); i++)
    for (unsigned int j = 0; j < rhs.ncols(); j++)
      tmp[i][j] = a - rhs[i][j];
	
  return tmp;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator-=(const qpMatrix<T>& rhs)
{
  if (m != rhs.ncols() || n != rhs.nrows())
    throw std::logic_error("Operator-=: matrices have different sizes");
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] -= rhs[i][j];
	
  return *this;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator-=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] -= a;
	
  return *this;
}

template <typename T>
qpMatrix<T> operator*(const qpMatrix<T>& lhs, const qpMatrix<T>& rhs)
{
  if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
    throw std::logic_error("Operator*: matrices have different sizes");
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] * rhs[i][j];
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator*(const qpMatrix<T>& lhs, const T& a)
{
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] * a;
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator*(const T& a, const qpMatrix<T>& rhs)
{
  qpMatrix<T> tmp(rhs.nrows(), rhs.ncols());
  for (unsigned int i = 0; i < rhs.nrows(); i++)
    for (unsigned int j = 0; j < rhs.ncols(); j++)
      tmp[i][j] = a * rhs[i][j];
	
  return tmp;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator*=(const qpMatrix<T>& rhs)
{
  if (m != rhs.ncols() || n != rhs.nrows())
    throw std::logic_error("Operator*=: matrices have different sizes");
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] *= rhs[i][j];
	
  return *this;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator*=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] *= a;
	
  return *this;
}

template <typename T>
qpMatrix<T> operator/(const qpMatrix<T>& lhs, const qpMatrix<T>& rhs)
{
  if (lhs.ncols() != rhs.ncols() || lhs.nrows() != rhs.nrows())
    throw std::logic_error("Operator+: matrices have different sizes");
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] / rhs[i][j];
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator/(const qpMatrix<T>& lhs, const T& a)
{
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = lhs[i][j] / a;
	
  return tmp;
}

template <typename T>
qpMatrix<T> operator/(const T& a, const qpMatrix<T>& rhs)
{
  qpMatrix<T> tmp(rhs.nrows(), rhs.ncols());
  for (unsigned int i = 0; i < rhs.nrows(); i++)
    for (unsigned int j = 0; j < rhs.ncols(); j++)
      tmp[i][j] = a / rhs[i][j];
	
  return tmp;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator/=(const qpMatrix<T>& rhs)
{
  if (m != rhs.ncols() || n != rhs.nrows())
    throw std::logic_error("Operator+=: matrices have different sizes");
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] /= rhs[i][j];
	
  return *this;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator/=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] /= a;
	
  return *this;
}

template <typename T>
qpMatrix<T> operator^(const qpMatrix<T>& lhs, const T& a)
{
  qpMatrix<T> tmp(lhs.nrows(), lhs.ncols());
  for (unsigned int i = 0; i < lhs.nrows(); i++)
    for (unsigned int j = 0; j < lhs.ncols(); j++)
      tmp[i][j] = pow(lhs[i][j], a);
	
  return tmp;
}

template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator^=(const qpMatrix<T>& rhs)
{
  if (m != rhs.ncols() || n != rhs.nrows())
    throw std::logic_error("Operator^=: matrices have different sizes");
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = pow(v[i][j], rhs[i][j]);
	
  return *this;
}


template <typename T>
inline qpMatrix<T>& qpMatrix<T>::operator^=(const T& a)
{
  for (unsigned int i = 0; i < n; i++)
    for (unsigned int j = 0; j < m; j++)
      v[i][j] = pow(v[i][j], a);
	
  return *this;
}

template <typename T>
inline qpMatrix<T>::operator qpVector<T>()
{
  if (n > 1 && m > 1)
    throw std::logic_error("Error matrix cast to vector: trying to cast a multi-dimensional matrix");
  if (n == 1)
    return extractRow(0);
  else
    return extractColumn(0);
}

template <typename T>
inline bool operator==(const qpMatrix<T>& a, const qpMatrix<T>& b)
{
  if (a.nrows() != b.nrows() || a.ncols() != b.ncols())
    throw std::logic_error("Matrices of different size are not confrontable");
  for (unsigned i = 0; i < a.nrows(); i++)
    for (unsigned j = 0; j < a.ncols(); j++)
      if (a[i][j] != b[i][j])
	return false;
  return true;
}

template <typename T>
inline bool operator!=(const qpMatrix<T>& a, const qpMatrix<T>& b)
{
  if (a.nrows() != b.nrows() || a.ncols() != b.ncols())
    throw std::logic_error("Matrices of different size are not confrontable");
  for (unsigned i = 0; i < a.nrows(); i++)
    for (unsigned j = 0; j < a.ncols(); j++)
      if (a[i][j] != b[i][j])
	return true;
  return false;
}



/**
   Input/Output 
*/
template <typename T>
std::ostream& operator<<(std::ostream& os, const qpMatrix<T>& m)
{
  os << std::endl << m.nrows() << " " << m.ncols() << std::endl;
  for (unsigned int i = 0; i < m.nrows(); i++)
    {
      for (unsigned int j = 0; j < m.ncols() - 1; j++)
	os << std::setw(20) << std::setprecision(16) << m[i][j] << ", ";
      os << std::setw(20) << std::setprecision(16) << m[i][m.ncols() - 1] << std::endl;
    }
	
  return os;
}

template <typename T>
std::istream& operator>>(std::istream& is, qpMatrix<T>& m)
{
  int rows, cols;
  char comma;
  is >> rows >> cols;
  m.resize(rows, cols);
  for (unsigned int i = 0; i < rows; i++)
    for (unsigned int j = 0; j < cols; j++)
      is >> m[i][j] >> comma;
	
  return is;
}

template <typename T>
T sign(const T& v)
{
  if (v >= (T)0.0)
    return (T)1.0;
  else
    return (T)-1.0;
}

template <typename T>
T dist(const T& a, const T& b)
{
  T abs_a = (T)fabs(a), abs_b = (T)fabs(b);
  if (abs_a > abs_b)
    return abs_a * sqrt((T)1.0 + (abs_b / abs_a) * (abs_b / abs_a));
  else
    return (abs_b == (T)0.0 ? (T)0.0 : abs_b * sqrt((T)1.0 + (abs_a / abs_b) * (abs_a / abs_b)));
}

template <typename T>
void svd(const qpMatrix<T>& A, qpMatrix<T>& U, qpVector<T>& W, qpMatrix<T>& V)
{
  int m = A.nrows(), n = A.ncols(), i, j, k, l, jj, nm;
  const unsigned int max_its = 30;
  bool flag;
  qpVector<T> rv1(n);
  U = A;
  W.resize(n);
  V.resize(n, n);
  T anorm, c, f, g, h, s, scale, x, y, z;
  g = scale = anorm = (T)0.0;
	
  // Householder reduction to bidiagonal form
  for (i = 0; i < n; i++)
    {
      l = i + 1;
      rv1[i] = scale * g;
      g = s = scale = (T)0.0;
      if (i < m)
	{
	  for (k = i; k < m; k++)
	    scale += fabs(U[k][i]);
	  if (scale != (T)0.0)
	    {
	      for (k = i; k < m; k++)
		{
		  U[k][i] /= scale;
		  s += U[k][i] * U[k][i];
		}
	      f = U[i][i];
	      g = -sign(f) * sqrt(s);
	      h = f * g - s;
	      U[i][i] = f - g;
	      for (j = l; j < n; j++)
		{
		  s = (T)0.0;
		  for (k = i; k < m; k++)
		    s += U[k][i] * U[k][j];
		  f = s / h;
		  for (k = i; k < m; k++)
		    U[k][j] += f * U[k][i];
		}
	      for (k = i; k < m; k++)
		U[k][i] *= scale;
	    }
	}
      W[i] = scale * g;
      g = s = scale = (T)0.0;
      if (i < m && i != n - 1)
	{
	  for (k = l; k < n; k++)
	    scale += fabs(U[i][k]);
	  if (scale != (T)0.0)
	    {
	      for (k = l; k < n; k++)
		{
		  U[i][k] /= scale;
		  s += U[i][k] * U[i][k];					
		}
	      f = U[i][l];
	      g = -sign(f) * sqrt(s);
	      h = f * g - s;
	      U[i][l] = f - g;
	      for (k = l; k <n; k++)
		rv1[k] = U[i][k] / h;
	      for (j = l; j < m; j++)
		{
		  s = (T)0.0;
		  for (k = l; k < n; k++)
		    s += U[j][k] * U[i][k];
		  for (k = l; k < n; k++)
		    U[j][k] += s * rv1[k];
		}
	      for (k = l; k < n; k++)
		U[i][k] *= scale;
	    }
	}
      anorm = ((std::max))(anorm, fabs(W[i]) + fabs(rv1[i]));
    }
  // Accumulation of right-hand transformations
  for (i = n - 1; i >= 0; i--)
    {
      if (i < n - 1) 
	{
	  if (g != (T)0.0)
	    {
	      for (j = l; j < n; j++)
		V[j][i] = (U[i][j] / U[i][l]) / g;
	      for (j = l; j < n; j++)
		{
		  s = (T)0.0;
		  for (k = l; k < n; k++)
		    s += U[i][k] * V[k][j];
		  for (k = l; k < n; k++)
		    V[k][j] += s * V[k][i];
		}
	    }
	  for (j = l; j < n; j++)
	    V[i][j] = V[j][i] = (T)0.0;
	}
      V[i][i] = (T)1.0;
      g = rv1[i];
      l = i;
    }
  // Accumulation of left-hand transformations
  for (i = ((std::min))(m, n) - 1; i >= 0; i--)
    {
      l = i + 1;
      g = W[i];
      for (j = l; j < n; j++)
	U[i][j] = (T)0.0;
      if (g != (T)0.0)
	{
	  g = (T)1.0 / g;
	  for (j = l; j < n; j++)
	    {
	      s = (T)0.0;
	      for (k = l; k < m; k++)
		s += U[k][i] * U[k][j];
	      f = (s / U[i][i]) * g;
	      for (k = i; k < m; k++)
		U[k][j] += f * U[k][i];
	    }
	  for (j = i; j < m; j++)
	    U[j][i] *= g;
	}
      else
	for (j = i; j < m; j++)
	  U[j][i] = (T)0.0;
      U[i][i]++;
    }
  // Diagonalization of the bidiagonal form: loop over singular values, and over allowed iterations.
  for (k = n - 1; k >= 0; k--)
    {
      for (unsigned int its = 0; its < max_its; its++)
	{
	  flag = true;
	  for (l = k; l >= 0; l--) // FIXME: in NR it was l >= 1 but there subscripts start from one
	    { // Test for splitting
	      nm = l - 1; // Note that rV[0] is always zero
	      if ((T)(fabs(rv1[l]) + anorm) == anorm)
		{
		  flag = false;
		  break;
		}
	      if ((T)(fabs(W[nm]) + anorm) == anorm)
		break;
	    }
	  if (flag)
	    {
	      // Cancellation of rv1[l], if l > 0 FIXME: it was l > 1 in NR
	      c = (T)0.0;
	      s = (T)1.0;
	      for (i = l; i <= k; i++)
		{
		  f = s * rv1[i];
		  rv1[i] *= c;
		  if ((T)(fabs(f) + anorm) == anorm)
		    break;
		  g = W[i];
		  h = dist(f, g);
		  W[i] = h;
		  h = (T)1.0 / h;
		  c = g * h;
		  s = -f * h;
		  for (j = 0; j < m; j++)
		    {
		      y = U[j][nm];
		      z = U[j][i];
		      U[j][nm] = y * c + z * s;
		      U[j][i] = z * c - y * s;
		    }
		}
	    }
	  z = W[k];
	  if (l == k)
	    {  // Convergence
	      if (z < (T)0.0)
		{ // Singular value is made nonnegative
		  W[k] = -z;
		  for (j = 0; j < n; j++)
		    V[j][k] = -V[j][k];
		}
	      break;
	    }
	  if (its == max_its)
	    throw std::logic_error("Error svd: no convergence in the maximum number of iterations");
	  x = W[l];
	  nm = k - 1;
	  y = W[nm];
	  g = rv1[nm];
	  h = rv1[k];
	  f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
	  g = dist(f, (T)1.0);
	  f = ((x - z) * (x + z) + h * ((y / (f + sign(f)*fabs(g))) - h)) / x;
	  c = s = (T)1.0; // Next QR transformation
	  for (j = l; j <= nm; j++)
	    {
	      i = j + 1;
	      g = rv1[i];
	      y = W[i];
	      h = s * g;
	      g *= c;
	      z = dist(f, h);
	      rv1[j] = z;
	      c = f / z;
	      s = h / z;
	      f = x * c + g * s;
	      g = g * c - x * s;
	      h = y * s;
	      y *= c;
	      for (jj = 0; jj < n; jj++)
		{
		  x = V[jj][j];
		  z = V[jj][i];
		  V[jj][j] = x * c + z * s;
		  V[jj][i] = z * c - x * s;
		}
	      z = dist(f, h);
	      W[j] = z; 
	      if (z != 0) // Rotation can be arbitrary if z = 0
		{
		  z = (T)1.0 / z;
		  c = f * z;
		  s = h * z;
		}
	      f = c * g + s * y;
	      x = c * y - s * g;
	      for (jj = 0; jj < m; jj++)
		{
		  y = U[jj][j];
		  z = U[jj][i];
		  U[jj][j] = y * c + z * s;
		  U[jj][i] = z * c - y * s;
		}
	    }
	  rv1[l] = (T)0.0;
	  rv1[k] = f;
	  W[k] = x;
	}
    }	
}

template <typename T>
qpMatrix<T> pinv(const qpMatrix<T>& A)
{
  qpMatrix<T> U, V, x, tmp(A.ncols(), A.nrows());
  qpVector<T> W;
  CanonicalBaseVector<T> e(0, A.nrows());
  svd(A, U, W, V);
  for (unsigned int i = 0; i < A.nrows(); i++)
    {
      e.reset(i);
      tmp.setColumn(i, dot_prod(dot_prod(dot_prod(V, qpMatrix<double>(DIAG, 1.0 / W, 0.0, W.size(), W.size())), t(U)), e));
    }
		
  return tmp;
}

template <typename T>
int lu(const qpMatrix<T>& A, qpMatrix<T>& LU, qpVector<unsigned int>& index)
{
  if (A.ncols() != A.nrows())
    throw std::logic_error("Error in LU decomposition: matrix must be squared");
  int i, p, j, k, n = A.ncols(), ex;
  T val, tmp;
  qpVector<T> d(n);
  LU = A;
  index.resize(n);
	
  ex = 1;
  for (i = 0; i < n; i++)
    {
      index[i] = i;
      val = (T)0.0;
      for (j = 0; j < n; j++)
	val = (std::max)(val, (T)fabs(LU[i][j]));
      if (val == (T)0.0)
	std::logic_error("Error in LU decomposition: matrix was singular");
      d[i] = val;
    }

  for (k = 0; k < n - 1; k++)
    {
      p = k;
      val = fabs(LU[k][k]) / d[k];
      for (i = k + 1; i < n; i++)
	{
	  tmp = fabs(LU[i][k]) / d[i];
	  if (tmp > val)
	    {
	      val = tmp;
	      p = i;
	    }
	}
      if (val == (T)0.0)
	std::logic_error("Error in LU decomposition: matrix was singular");
      if (p > k)
	{
	  ex = -ex;
	  std::swap(index[k], index[p]);
	  std::swap(d[k], d[p]);
	  for (j = 0; j < n; j++)
	    std::swap(LU[k][j], LU[p][j]);
	}
		
      for (i = k + 1; i < n; i++)
	{
	  LU[i][k] /= LU[k][k];
	  for (j = k + 1; j < n; j++)
	    LU[i][j] -= LU[i][k] * LU[k][j];
	}
    }
  if (LU[n - 1][n - 1] == (T)0.0)
    std::logic_error("Error in LU decomposition: matrix was singular");
		
  return ex;
}

template <typename T>
qpVector<T> lu_solve(const qpMatrix<T>& LU, const qpVector<T>& b, qpVector<unsigned int>& index)
{
  if (LU.ncols() != LU.nrows())
    throw std::logic_error("Error in LU solve: LU matrix should be squared");
  unsigned int n = LU.ncols();
  if (b.size() != n)
    throw std::logic_error("Error in LU solve: b vector must be of the same dimensions of LU matrix");
  qpVector<T> x((T)0.0, n);
  int i, j, p;
  T sum;
	
  p = index[0];
  x[0] = b[p];
	
  for (i = 1; i < n; i++)
    {
      sum = (T)0.0;
      for (j = 0; j < i; j++)
	sum += LU[i][j] * x[j];
      p = index[i];
      x[i] = b[p] - sum;
    }
  x[n - 1] /= LU[n - 1][n - 1];
  for (i = n - 2; i >= 0; i--)
    {
      sum = (T)0.0;
      for (j = i + 1; j < n; j++)
	sum += LU[i][j] * x[j];
      x[i] = (x[i] - sum) / LU[i][i];
    }
  return x;
}

template <typename T>
void lu_solve(const qpMatrix<T>& LU, qpVector<T>& x, const qpVector<T>& b, qpVector<unsigned int>& index)
{
  x = lu_solve(LU, b, index);
}

template <typename T>
qpMatrix<T> lu_inverse(const qpMatrix<T>& A)
{
  if (A.ncols() != A.nrows())
    throw std::logic_error("Error in LU invert: matrix must be squared");	
  unsigned int n = A.ncols();
  qpMatrix<T> A1(n, n), LU;
  qpVector<unsigned int> index;
	
  lu(A, LU, index);
  CanonicalBaseVector<T> e(0, n);
  for (unsigned i = 0; i < n; i++)
    {
      e.reset(i);
      A1.setColumn(i, lu_solve(LU, e, index));
    }
	
  return A1;
}

template <typename T>
T lu_det(const qpMatrix<T>& A)
{
  if (A.ncols() != A.nrows())
    throw std::logic_error("Error in LU determinant: matrix must be squared");	
  unsigned int d;
  qpMatrix<T> LU;
  qpVector<unsigned int> index;
	
  d = lu(A, LU, index);
	
  return d * prod(LU.extractDiag());
}

template <typename T>
void cholesky(const qpMatrix<T> A, qpMatrix<T>& LL) 
{
  if (A.ncols() != A.nrows())
    throw std::logic_error("Error in Cholesky decomposition: matrix must be squared");
  int n = A.ncols();
  double sum;
  LL = A;
	
  for (unsigned int i = 0; i < n; i++)
    {
      for (unsigned int j = i; j < n; j++)
	{
	  sum = LL[i][j];
	  for (int k = i - 1; k >= 0; k--)
	    sum -= LL[i][k] * LL[j][k];
	  if (i == j) 
	    {
	      if (sum <= 0.0)
		throw std::logic_error("Error in Cholesky decomposition: matrix is not postive definite");
	      LL[i][i] = sqrt(sum);
	    }
	  else
	    LL[j][i] = sum / LL[i][i];
	}
      for (unsigned int k = i + 1; k < n; k++)
	LL[i][k] = LL[k][i];
    } 
}

template <typename T>
qpMatrix<T> cholesky(const qpMatrix<T> A) 
{
  qpMatrix<T> LL;
  cholesky(A, LL);
	
  return LL;
}

template <typename T>
qpVector<T> cholesky_solve(const qpMatrix<T>& LL, const qpVector<T>& b)
{
  if (LL.ncols() != LL.nrows())
    throw std::logic_error("Error in Cholesky solve: matrix must be squared");
  unsigned int n = LL.ncols();
  if (b.size() != n)
    throw std::logic_error("Error in Cholesky decomposition: b vector must be of the same dimensions of LU matrix");
  qpVector<T> x, y;
	
  /* Solve L * y = b */
  forward_elimination(LL, y, b);
  /* Solve L^T * x = y */
  backward_elimination(LL, x, y);
	
  return x;
}

template <typename T>
void cholesky_solve(const qpMatrix<T>& LL, qpVector<T>& x, const qpVector<T>& b)
{
  x = cholesky_solve(LL, b);
}

template <typename T>
void forward_elimination(const qpMatrix<T>& L, qpVector<T>& y, const qpVector<T> b)
{
  if (L.ncols() != L.nrows())
    throw std::logic_error("Error in Forward elimination: matrix must be squared (lower triangular)");
  if (b.size() != L.nrows())
    throw std::logic_error("Error in Forward elimination: b vector must be of the same dimensions of L matrix");
  unsigned int n = b.size();
  y.resize(n);
	
  y[0] = b[0] / L[0][0];
  for (unsigned int i = 1; i < n; i++)
    {
      y[i] = b[i];
      for (unsigned int j = 0; j < i; j++)
	y[i] -= L[i][j] * y[j];
      y[i] = y[i] / L[i][i];
    }
}

template <typename T>
qpVector<T> forward_elimination(const qpMatrix<T>& L, const qpVector<T> b)
{
  qpVector<T> y;
  forward_elimination(L, y, b);
	
  return y;
}

template <typename T>
void backward_elimination(const qpMatrix<T>& U, qpVector<T>& x, const qpVector<T>& y)
{
  if (U.ncols() != U.nrows())
    throw std::logic_error("Error in Backward elimination: matrix must be squared (upper triangular)");
  if (y.size() != U.nrows())
    throw std::logic_error("Error in Backward elimination: b vector must be of the same dimensions of U matrix");
  int n = y.size();
  x.resize(n);
	
  x[n - 1] = y[n - 1] / U[n - 1][n - 1];
  for (int i = n - 2; i >= 0; i--)
    {
      x[i] = y[i];
      for (int j = i + 1; j < n; j++)
	x[i] -= U[i][j] * x[j];
      x[i] = x[i] / U[i][i];
    }
}

template <typename T>
qpVector<T> backward_elimination(const qpMatrix<T>& U, const qpVector<T> y)
{
  qpVector<T> x;
  forward_elimination(U, x, y);
	
  return x;
}

/* Setting default linear systems machinery */

#define det lu_det
#define inverse lu_inverse
#define solve lu_solve

/* Random */

template <typename T>
void random(qpMatrix<T>& m)
{
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      m[i][j] = (T)(rand() / double(RAND_MAX));
}

/**
   Aggregate functions
*/

template <typename T>
qpVector<T> sum(const qpMatrix<T>& m)
{
  qpVector<T> tmp((T)0, m.ncols());
  for (unsigned int j = 0; j < m.ncols(); j++)
    for (unsigned int i = 0; i < m.nrows(); i++)
      tmp[j] += m[i][j];
  return tmp;
}

template <typename T>
qpVector<T> r_sum(const qpMatrix<T>& m)
{
  qpVector<T> tmp((T)0, m.nrows());
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp[i] += m[i][j];
  return tmp;
}

template <typename T>
T all_sum(const qpMatrix<T>& m)
{
  T tmp = (T)0;
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp += m[i][j];
  return tmp;
}

template <typename T>
qpVector<T> prod(const qpMatrix<T>& m)
{
  qpVector<T> tmp((T)1, m.ncols());
  for (unsigned int j = 0; j < m.ncols(); j++)
    for (unsigned int i = 0; i < m.nrows(); i++)
      tmp[j] *= m[i][j];
  return tmp;
}

template <typename T>
qpVector<T> r_prod(const qpMatrix<T>& m)
{
  qpVector<T> tmp((T)1, m.nrows());
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp[i] *= m[i][j];
  return tmp;
}

template <typename T>
T all_prod(const qpMatrix<T>& m)
{
  T tmp = (T)1;
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp *= m[i][j];
  return tmp;
}

template <typename T>
qpVector<T> mean(const qpMatrix<T>& m)
{
  qpVector<T> res((T)0, m.ncols());
  for (unsigned int j = 0; j < m.ncols(); j++)
    {
      for (unsigned int i = 0; i < m.nrows(); i++)
	res[j] += m[i][j];
      res[j] /= m.nrows();
    }
	
  return res;
}

template <typename T>
qpVector<T> r_mean(const qpMatrix<T>& m)
{
  qpVector<T> res((T)0, m.rows());
  for (unsigned int i = 0; i < m.nrows(); i++)
    {
      for (unsigned int j = 0; j < m.ncols(); j++)
	res[i] += m[i][j];
      res[i] /= m.nrows();
    }
	
  return res;
}

template <typename T>
T all_mean(const qpMatrix<T>& m)
{
  T tmp = (T)0;
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp += m[i][j];
  return tmp / (m.nrows() * m.ncols());
}

template <typename T>
qpVector<T> var(const qpMatrix<T>& m, bool sample_correction = false)
{
  qpVector<T> res((T)0, m.ncols());
  unsigned int n = m.nrows();
  double sum, ssum;
  for (unsigned int j = 0; j < m.ncols(); j++)
    {	
      sum = (T)0.0; ssum = (T)0.0;
      for (unsigned int i = 0; i < m.nrows(); i++)
	{
	  sum += m[i][j];
	  ssum += (m[i][j] * m[i][j]);
	}
      if (!sample_correction)
	res[j] = (ssum / n) - (sum / n) * (sum / n);
      else
	res[j] = n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);		 
    }
	
  return res;
}

template <typename T>
qpVector<T> stdev(const qpMatrix<T>& m, bool sample_correction = false)
{
  return vec_sqrt(var(m, sample_correction));
}

template <typename T>
qpVector<T> r_var(const qpMatrix<T>& m, bool sample_correction = false)
{
  qpVector<T> res((T)0, m.nrows());
  double sum, ssum;
  unsigned int n = m.ncols();
  for (unsigned int i = 0; i < m.nrows(); i++)
    {	
      sum = 0.0; ssum = 0.0;
      for (unsigned int j = 0; j < m.ncols(); j++)
	{
	  sum += m[i][j];
	  ssum += (m[i][j] * m[i][j]);
	}
      if (!sample_correction)
	res[i] = (ssum / n) - (sum / n) * (sum / n);
      else
	res[i] = n * ((ssum / n) - (sum / n) * (sum / n)) / (n - 1);
    }
	
  return res;
}

template <typename T>
qpVector<T> r_stdev(const qpMatrix<T>& m, bool sample_correction = false)
{
  return vec_sqrt(r_var(m, sample_correction));
}

//template <typename T>
//qpVector<T> max(const qpMatrix<T>& m)
//{
//  qpVector<T> res(m.ncols());
//  double value;
//  for (unsigned int j = 0; j < m.ncols(); j++)
//    {
//      value = m[0][j];
//      for (unsigned int i = 1; i < m.nrows(); i++)
//	value = (std::max)(m[i][j], value);
//      res[j] = value;
//    }
//	
//  return res;
//}

template <typename T>
qpVector<T> r_max(const qpMatrix<T>& m)
{
  qpVector<T> res(m.nrows());
  double value;
  for (unsigned int i = 0; i < m.nrows(); i++)
    {
      value = m[i][0];
      for (unsigned int j = 1; j < m.ncols(); j++)
	value = (std::max)(m[i][j], value);
      res[i] = value;
    }
	
  return res;
}

//template <typename T>
//qpVector<T> min(const qpMatrix<T>& m)
//{
//  qpVector<T> res(m.ncols());
//  double value;
//  for (unsigned int j = 0; j < m.ncols(); j++)
//    {
//      value = m[0][j];
//      for (unsigned int i = 1; i < m.nrows(); i++)
//	value = (std::min)(m[i][j], value);
//      res[j] = value;
//    }
//	
//  return res;
//}

template <typename T>
qpVector<T> r_min(const qpMatrix<T>& m)
{
  qpVector<T> res(m.nrows());
  double value;
  for (unsigned int i = 0; i < m.nrows(); i++)
    {
      value = m[i][0];
      for (unsigned int j = 1; j < m.ncols(); j++)
	value = (std::min)(m[i][j], value);
      res[i] = value;
    }
	
  return res;
}



/**
   Single element mathematical functions
*/

template <typename T>
qpMatrix<T> exp(const qpMatrix<T>&m)
{
  qpMatrix<T> tmp(m.nrows(), m.ncols());
	
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp[i][j] = exp(m[i][j]);
	
  return tmp;
}

template <typename T>
qpMatrix<T> mat_sqrt(const qpMatrix<T>&m)
{
  qpMatrix<T> tmp(m.nrows(), m.ncols());
	
  for (unsigned int i = 0; i < m.nrows(); i++)
    for (unsigned int j = 0; j < m.ncols(); j++)
      tmp[i][j] = sqrt(m[i][j]);
	
  return tmp;
}

/**
   qpMatrix operators
*/

template <typename T>
qpMatrix<T> kron(const qpVector<T>& b, const qpVector<T>& a)
{
  qpMatrix<T> tmp(b.size(), a.size());
  for (unsigned int i = 0; i < b.size(); i++)
    for (unsigned int j = 0; j < a.size(); j++)
      tmp[i][j] = a[j] * b[i];
	
  return tmp;
}

template <typename T>
qpMatrix<T> t(const qpMatrix<T>& a)
{
  qpMatrix<T> tmp(a.ncols(), a.nrows());
  for (unsigned int i = 0; i < a.nrows(); i++)
    for (unsigned int j = 0; j < a.ncols(); j++)
      tmp[j][i] = a[i][j];
	
  return tmp;
}

template <typename T>
qpMatrix<T> dot_prod(const qpMatrix<T>& a, const qpMatrix<T>& b)
{
  if (a.ncols() != b.nrows())
    throw std::logic_error("Error matrix dot product: dimensions of the matrices are not compatible");
  qpMatrix<T> tmp(a.nrows(), b.ncols());
  for (unsigned int i = 0; i < tmp.nrows(); i++)
    for (unsigned int j = 0; j < tmp.ncols(); j++)
      {
	tmp[i][j] = (T)0;
	for (unsigned int k = 0; k < a.ncols(); k++)
	  tmp[i][j] += a[i][k] * b[k][j];
      }
			
  return tmp;
}

template <typename T>
qpMatrix<T> dot_prod(const qpMatrix<T>& a, const qpVector<T>& b)
{
  if (a.ncols() != b.size())
    throw std::logic_error("Error matrix dot product: dimensions of the matrix and the vector are not compatible");
  qpMatrix<T> tmp(a.nrows(), 1);
  for (unsigned int i = 0; i < tmp.nrows(); i++)
    {
      tmp[i][0] = (T)0;
      for (unsigned int k = 0; k < a.ncols(); k++)
	tmp[i][0] += a[i][k] * b[k];
    }
		
  return tmp;
}

template <typename T>
qpMatrix<T> dot_prod(const qpVector<T>& a, const qpMatrix<T>& b)
{
  if (a.size() != b.nrows())
    throw std::logic_error("Error matrix dot product: dimensions of the vector and matrix are not compatible");
  qpMatrix<T> tmp(1, b.ncols());
  for (unsigned int j = 0; j < tmp.ncols(); j++)
    {
      tmp[0][j] = (T)0;
      for (unsigned int k = 0; k < a.size(); k++)
	tmp[0][j] += a[k] * b[k][j];
    }
		
  return tmp;
}

template <typename T>
inline qpMatrix<double> rank(const qpMatrix<T> m)
{
  qpMatrix<double> tmp(m.nrows(), m.ncols());
  for (unsigned int j = 0; j < m.ncols(); j++)
    tmp.setColumn(j, rank<T>(m.extractColumn(j)));
  
  return tmp;                  
}

template <typename T>
inline qpMatrix<double> r_rank(const qpMatrix<T> m)
{
  qpMatrix<double> tmp(m.nrows(), m.ncols());
  for (unsigned int i = 0; i < m.nrows(); i++)
    tmp.setRow(i, rank<T>(m.extractRow(i)));
  
  return tmp;                  
}

} // namespace quadprogpp

#endif // define _ARRAY_HH_
