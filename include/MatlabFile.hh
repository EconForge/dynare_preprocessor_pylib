/*
 * Copyright (C) 2009 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */
/*
  Usefull documentation: Matlab 7 Mat-File Format
  -----------------------------------------------
  revision: October 2008 PDF only Rereleased for Version 7.7 (Release 2008b)
  available at: http://www.mathworks.com/access/helpdesk/help/pdf_doc/matlab/matfile_format.pdf
*/

#ifndef _MAT_FILE_HH
#define _MAT_FILE_HH

#include <cstdio>
#include <cstdlib>
#include <string>
#include <fstream>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
#include <math.h>
//! zlib needed to uncompress the mat-file. It is available with GCC 4.3.2 but it needs a dll !!
//! => to avoid compress MatFile, save is used with option '-v6' in save_params_and_strady_state.m
//#include "zlib.h"
using namespace std;


enum Data_Type
  {
    miINT8       = 1,  //8 bit, signed
    miUINT8      = 2,  //8 bit, unsigned
    miINT16      = 3,  //16-bit, signed
    miUINT16     = 4,  //16-bit, unsigned
    miINT32      = 5,  //32-bit, signed
    miUINT32     = 6,  //32-bit, unsigned
    miSINGLE     = 7,  //IEEE® 754 single format
    miDOUBLE     = 9,  //IEEE 754 double format
    miINT64      = 12, //64-bit, signed
    miUINT64     = 13, //64-bit, unsigned
    miMATRIX     = 14, //MATLAB array
    miCOMPRESSED = 15, //Compressed Data
    miUTF8       = 16, //Unicode UTF-8 Encoded Character Data
    miUTF16      = 17, //Unicode UTF-16 Encoded Character Data
    miUTF32      = 18  //Unicode UTF-32 Encoded Character Data
  };

enum Array_Type
  {
    Cell_array              = 1,
    Structure_              = 2,
    Object_                 = 3,
    Character_array         = 4,
    Sparse_array            = 5,
    Double_precision_array  = 6,
    Single_precision_array  = 7,
    Signed_integer_8_bit    = 8,
    Unsigned_integer_8_bit  = 9,
    Signed_integer_16_bit   = 10,
    Unsigned_integer_16_bit = 11,
    Signed_integer_32_bit   = 12,
    Unsigned_integer_32_bit = 13
  };


enum Function_Returned_Type
  {
    Numerical    =1,
    AlphaNumeric =2,
    Matrix       =3,
    Compressed   =4,
    Unknown      =5
  };

class ArrayElem;
class SimpleElem;

typedef long long LongLongInt;
typedef long int LongInt;
typedef short int Int;
//typedef char ShortInt;
typedef unsigned long int uLongInt ;
typedef unsigned short int uShortInt ;
typedef short int ShortInt ;
typedef class SimpleElem *PSimpleElem;

//!Header of MatFile
typedef struct Header
{
  char Theader[124];
  short int Version;
  char Edian_Indicator[2];
}
  Header_t;

typedef struct Data_Header
{
  ShortInt S_Number_of_Bytes;
  ShortInt DataType;
  uLongInt Number_of_Bytes;
}
  Data_Header_t;

typedef struct Array_Flag
{
  Data_Header_t tag;
  unsigned char classe;
  unsigned char flag;
  char undef1[2];
  uLongInt nzmax;
}
  Array_Flag_t;


typedef struct returned_ReadData
{
  SimpleElem* Simple;
  Function_Returned_Type Type;
} returned_ReadData_t;


typedef struct FlagStructure
{
  bool no_name;
  bool character;
} FlagStructure_t;

typedef struct CollectStruct
{
  /*vector<string> variable_name;
    vector< vector<double> > variable_double;
    vector< vector<string> > variable_string;*/
  string tmp_name;
  map<string,vector<string> > variable_string_name;
  map<string,vector<double> > variable_double_name;
}
  CollectStruct;

typedef vector<int> Array_Dimensions_t;





//! Base class for simple elements in Mat-File
class SimpleElem
{
public:
  bool verbose;
  vector<double> VNumeric;
  vector<string> Vstr;
  ArrayElem *array_elem;
  int Type;
  SimpleElem();
  virtual ~SimpleElem();
  virtual double ReadNum(char* InBuff, int* pBuff) const{return(0);};
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const{return(NULL);};
  virtual Data_Header_t ReadDataHeader(char* InBuff, int* pBuff) const;
  virtual int size() const{cout << "oups\n";return(0);};
  void Print() const;
  void Delete() const;
  void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  returned_ReadData_t ReadData(char* InBuff, int* pBuff) const;
  returned_ReadData_t Get_Data_Class(Data_Header_t data_header) const;
  void DataProceed(Data_Header data_header, char* InBuff, int* pBuff, FlagStructure flag);
};


class UTF8 : public SimpleElem
{
public:
  virtual int size() const {return(1);};
  virtual double ReadNum(char* InBuff, int* pBuff) const {return(NULL);};
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const;
};

class UTF16 : public SimpleElem
{
public:
  virtual int size() const {return(2);};
  virtual double ReadNum(char* InBuff, int* pBuff) const {return(NULL);};
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const;
};

class UTF32 : public SimpleElem
{
public:
  virtual int size() const {return(4);};
  virtual double ReadNum(char* InBuff, int* pBuff) const {return(NULL);};
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const;
};

class INT8 : public SimpleElem
{
public:
  virtual int size() const {return(1);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};

class INT16 : public SimpleElem
{
public:
  virtual int size() const {return(2);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};


class INT32 : public SimpleElem
{
public:
  virtual int size() const {return(4);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};


class INT64 : public SimpleElem
{
public:
  virtual int size() const {return(8);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};


class Single : public SimpleElem
{
public:
  virtual int size() const {return(4);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};


class Double : public SimpleElem
{
public:
  virtual int size() const {return(8);};
  virtual double ReadNum(char* InBuff, int* pBuff) const;
  virtual string ReadAlph(char* InBuff, int* pBuff, int Size) const {return(NULL);};
};

//! Base class for Array Element in Mat-File
class ArrayElem : public SimpleElem
{
private:

protected:

public:
  int Cell_number, Structure_number, Matrix_Elem_number;
  Array_Type Type;
  vector<PSimpleElem> VCell;
  vector<string> Structure_Elem_name;
  bool array_complex, array_global, array_logical;
  string variable_name;
  int array_nzmax;
  int number_of_dimensions;
  vector<int> dimension;
  vector<double> Double_value;
  vector<string> String_value;

  Array_Type type;
  ArrayElem();
  virtual LongInt ReadINT32(char* InBuff, int* pBuff) const;
  virtual Array_Flag_t ReadArrayFlag(char* InBuff, int* pBuff) /*const*/;
  virtual void ReadArrayDimension(char* InBuff, int* pBuff) /*const*/;
  virtual void ReadArrayName(char* InBuff, int* pBuff) /*const*/;
  virtual void ReadStructureNames(char* InBuff, int* pBuff);
  virtual ArrayElem* ReadArray_class(char* InBuff, int* pBuff) const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag) /*const*/ {cout << "oups..\n";};
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual ~ArrayElem();
};

class CellArray : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Structure : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Object : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class CharacterArray : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class SparseArray : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class DoublePrecisionArray : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class SinglePrecisionArray : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit8SignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit8UnsignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit16SignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit16UnsignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit32SignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};

class Bit32UnsignedInteger : public ArrayElem
{
public:
  virtual void Collect(const string &name, bool found, CollectStruct &collect_struct) const;
  virtual void Print() const;
  virtual void Delete() const;
  virtual void ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag);
};


class MatlabFile
{
public:
  Header_t header;
  vector<PSimpleElem> VSimpl;
  MatlabFile();
  ~MatlabFile();
  void MatFileRead(string filename);
  void MatFilePrint();
  void Delete();
  bool Collect(const string &name, CollectStruct &collect_struct) const;
  Data_Header_t ReadDataHeader(ifstream &MatFile);
};
#endif
