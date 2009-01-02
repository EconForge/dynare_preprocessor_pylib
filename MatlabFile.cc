/*
 * Copyright (C) 2006-2008 Dynare Team
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
#include "MatlabFile.hh"


MatlabFile::MatlabFile()
{
}

MatlabFile::~MatlabFile()
{
}

ArrayElem::ArrayElem()
{
}

ArrayElem::~ArrayElem()
{
}


string
UTF8::ReadAlph(char* InBuff, int* pBuff, int Size) const
{
  char tmp_c[Size+1];
  string str;
  memcpy(&tmp_c, InBuff+*pBuff, Size);
  tmp_c[Size]=0;
  str.assign(tmp_c);
  if(Size<=4)
    *pBuff += 4;
  else if(Size % 8)
    *pBuff += 8*ceil(double(Size) / 8);
  else
    *pBuff += Size;
  return(str);
}


string
UTF16::ReadAlph(char* InBuff, int* pBuff, int Size) const
{
  string str("");
  for(int i=0;i<Size;i++)
    str.append(&(InBuff[i*2+*pBuff]));
  if(Size*2<=4)
    *pBuff += 4;
  else if((Size*2) % 8)
    *pBuff += 8*ceil(double(Size*2) / 8);
  else
    *pBuff += Size*2;
  return(str);
}

string
UTF32::ReadAlph(char* InBuff, int* pBuff, int Size) const
{
  string str("");
  for(int i=0;i<Size;i++)
    str.append(&(InBuff[i*4+*pBuff]));
  if((Size*4) % 8)
    *pBuff += 8*ceil(double(Size*4) / 8);
  else
    *pBuff += Size*4;
  return(str);
}

double
INT8::ReadNum(char* InBuff, int* pBuff) const
{
  char val;
  val = InBuff[*pBuff];
  *pBuff += sizeof(val);
  return(val);
}


double
INT16::ReadNum(char* InBuff, int* pBuff) const
{
  Int val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}

double
INT32::ReadNum(char* InBuff, int* pBuff) const
{
  LongInt val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}


double
INT64::ReadNum(char* InBuff, int* pBuff) const
{
  LongLongInt val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}

double
Single::ReadNum(char* InBuff, int* pBuff) const
{
  float val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}

double
Double::ReadNum(char* InBuff, int* pBuff) const
{
  double val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}



Data_Header_t
SimpleElem::ReadDataHeader(char* InBuff, int* pBuff) const
{
  Data_Header_t data_header;
  memcpy(&data_header.DataType, InBuff+*pBuff, sizeof(data_header.DataType));
  *pBuff += 2;
  memcpy(&data_header.S_Number_of_Bytes, InBuff+*pBuff, sizeof(data_header.S_Number_of_Bytes));
  *pBuff += 2;
  if(data_header.S_Number_of_Bytes!=0)
    data_header.Number_of_Bytes=data_header.S_Number_of_Bytes;
  else
    {
      memcpy(&data_header.Number_of_Bytes, InBuff+*pBuff, sizeof(data_header.Number_of_Bytes));
      *pBuff += 4;
    }
  return(data_header);
}


SimpleElem::SimpleElem()
{
  verbose=true;
  array_elem=NULL;
  Type=-1;
}

SimpleElem::~SimpleElem()
{
}

returned_ReadData_t
SimpleElem::Get_Data_Class(Data_Header_t data_header) const
{
  returned_ReadData_t ret;
  switch(data_header.DataType)
    {
      case   miINT8:       //8 bit, signed
      case   miUINT8:      //8 bit, unsigned
        ret.Type = Numerical;
        ret.Simple = new INT8;
        return(ret);
        break;
      case   miINT16:      //16-bit, signed
      case   miUINT16:     //16-bit, unsigned
        ret.Type = Numerical;
        ret.Simple = new INT16;
        return(ret);
        break;
      case   miINT32:      //32-bit, signed
      case   miUINT32:     //32-bit, unsigned
        ret.Type = Numerical;
        ret.Simple = new INT32;
        return(ret);
        break;
      case   miINT64:      //64-bit, signed
      case   miUINT64:     //64-bit, unsigned
        ret.Type = Numerical;
        ret.Simple = new INT64;
        return(ret);
        break;
      case   miSINGLE:     //IEEE® 754 single format
        ret.Type = Numerical;
        ret.Simple = new Single;
        return(ret);
        break;
      case   miDOUBLE:     //IEEE 754 double format
        ret.Type = Numerical;
        ret.Simple = new Double;
        return(ret);
        break;
      case   miMATRIX:     //MATLAB array
        ret.Type = Matrix;
        ret.Simple = NULL;
        return(ret);
        break;
      case   miCOMPRESSED: //Compressed Data
        ret.Type = Compressed;
        ret.Simple = NULL;
        return(ret);
        break;
      case   miUTF8:       //Unicode UTF-8 Encoded Character Data
        ret.Type = AlphaNumeric;
        ret.Simple = new UTF8;
        return(ret);
        break;
      case   miUTF16:      //Unicode UTF-16 Encoded Character Data
        ret.Type = AlphaNumeric;
        ret.Simple = new UTF16;
        return(ret);
        break;
      case   miUTF32:      //Unicode UTF-32 Encoded Character Data
        ret.Type = AlphaNumeric;
        ret.Simple = new UTF32;
        return(ret);
        break;
      default:
        ret.Type = Unknown;
        ret.Simple = NULL;
        return(ret);
    }
}

void
SimpleElem::DataProceed(Data_Header data_header, char* InBuff, int* pBuff, FlagStructure_t flag)
{
  ArrayElem matrix;
  returned_ReadData ret;
  ret = Get_Data_Class(data_header);
  Type = ret.Type;
  double tmpv;
  switch(ret.Type)
    {
      case Numerical:
        if(data_header.Number_of_Bytes/ret.Simple->size())
          for(unsigned int i=0;i<data_header.Number_of_Bytes/ret.Simple->size();i++)
            {
              tmpv=ret.Simple->ReadNum(InBuff, pBuff);
              VNumeric.push_back(tmpv);
            }
        //to align pBuff on a 4 Bytes
        if(*pBuff % 4)
          *pBuff += 4-(*pBuff % 4);
        delete ret.Simple;
        break;
      case AlphaNumeric:
        for(unsigned int i=0;i<data_header.Number_of_Bytes/ret.Simple->size();i++)
          Vstr.push_back(ret.Simple->ReadAlph(InBuff, pBuff, data_header.Number_of_Bytes));
        //to align pBuff on a 4 Bytes
        if(*pBuff % 4)
          *pBuff += 4-(*pBuff % 4);
        delete ret.Simple;
        break;
      case Matrix:
        array_elem = matrix.ReadArray_class(InBuff, pBuff);
        array_elem->ReadArray(InBuff, pBuff, flag);
        break;
      case Compressed:
        cerr << "Error: Compressed data in Mat-file not implemnted yet!\n set option -v6 when saving to Mat-file\n";
        exit(EXIT_FAILURE);
      case Unknown:
        cerr << "Error: Mat-file format use incomptible format with Matlab 7 specification!\n set option -v6 when saving to Mat-file\n";
        exit(EXIT_FAILURE);
    }
}

returned_ReadData_t
SimpleElem::ReadData(char* InBuff, int* pBuff) const
{
  Data_Header_t data_header;
  data_header = ReadDataHeader(InBuff, pBuff);
  return(Get_Data_Class(data_header));
}

void
SimpleElem::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
  if(VNumeric.size())
    {
      if(found)
        collect_struct.variable_double_name[collect_struct.tmp_name]=VNumeric;
    }
  else if(Vstr.size())
    {
      if(found)
        collect_struct.variable_string_name[collect_struct.tmp_name]=Vstr;
    }
  else
    {
      if(array_elem)
        {
          array_elem->Collect(name, found, collect_struct);
        }
   }
}


void
SimpleElem::Print() const
{
  if(VNumeric.size())
    {
      for(vector<double>::const_iterator it=VNumeric.begin(); it!=VNumeric.end(); it++)
        cout << " " << *it;
    }
  else if(Vstr.size())
    {
      for(vector<string>::const_iterator it=Vstr.begin(); it!=Vstr.end(); it++)
        cout << " " << *it;
    }
  else
    {
      if(array_elem)
        array_elem->Print();
    }
}

void
SimpleElem::Delete() const
{
  if(array_elem)
    {
      array_elem->Delete();
      delete array_elem;
    }
}

LongInt
ArrayElem::ReadINT32(char* InBuff, int* pBuff) const
{
  LongInt val;
  memcpy(&val, InBuff+*pBuff, sizeof(val));
  *pBuff += sizeof(val);
  return(val);
}


Array_Flag_t
ArrayElem::ReadArrayFlag(char* InBuff, int* pBuff) /*const*/
{
  Array_Flag_t array_flag;
  memcpy(&array_flag, InBuff+*pBuff, sizeof(array_flag));
  *pBuff += sizeof(array_flag);
  array_complex = (bool)(array_flag.flag & 16);
  array_global = array_flag.flag & 32;
  array_logical = array_flag.flag & 64;
  type = (Array_Type)array_flag.classe;
  if(type==Sparse_array)
    array_nzmax = array_flag.nzmax;
  return(array_flag);
}

void
ArrayElem::ReadArrayDimension(char* InBuff, int* pBuff) /*const*/
{
  Data_Header_t data_header;
  data_header = ReadDataHeader(InBuff, pBuff);
  number_of_dimensions = data_header.Number_of_Bytes/4;
  for(int i=0; i<number_of_dimensions; i++)
    {
      double tmp_d=ReadINT32(InBuff, pBuff);
      dimension.push_back(tmp_d);
    }
  if(number_of_dimensions % 2)
    *pBuff += sizeof(dimension[0]);
}

void
ArrayElem::ReadStructureNames(char* InBuff, int* pBuff)
{
  Data_Header_t data_header, data_header_2;
  LongInt Field_name_length;
  data_header=ReadDataHeader(InBuff, pBuff);
  Field_name_length=ReadINT32(InBuff, pBuff);
  data_header_2=ReadDataHeader(InBuff, pBuff);
  Structure_number = data_header_2.Number_of_Bytes/Field_name_length;
  char tmp_c[Field_name_length];
  for(int i=0; i<Structure_number;i++)
    {
      memcpy(tmp_c, InBuff+*pBuff, Field_name_length);
      *pBuff += Field_name_length;
      string variable_name(tmp_c);
      Structure_Elem_name.push_back(variable_name);
    }
}

void
ArrayElem::ReadArrayName(char* InBuff, int* pBuff) /*const*/
{
  Data_Header_t data_header;
  data_header = ReadDataHeader(InBuff, pBuff);
  char tmp_c[data_header.Number_of_Bytes+1];
  memcpy(&tmp_c, InBuff+*pBuff, data_header.Number_of_Bytes);
  tmp_c[data_header.Number_of_Bytes]=0;
  variable_name.assign(tmp_c);
  if(data_header.Number_of_Bytes<=4)
    *pBuff += 4;
  else if(data_header.Number_of_Bytes % 8)
    *pBuff += 8*ceil(double(data_header.Number_of_Bytes) / 8);
  else
    *pBuff += data_header.Number_of_Bytes;
}


ArrayElem*
ArrayElem::ReadArray_class(char* InBuff, int* pBuff) const
{
  Array_Flag_t array_flag;
  ArrayElem array_elem;
  array_flag=array_elem.ReadArrayFlag(InBuff, pBuff);
  switch(array_flag.classe)
    {
       case Cell_array:
         return(new CellArray);
         break;
       case Structure_:
         return(new Structure);
         break;
       case Object_:
         return(new Object);
         break;
       case Character_array:
         return(new CharacterArray);
         break;
       case Sparse_array:
         return(new SparseArray);
         break;
       case Double_precision_array:
         return(new DoublePrecisionArray);
         break;
       case Single_precision_array:
         return(new SinglePrecisionArray);
         break;
       case Signed_integer_8_bit:
         return(new Bit8SignedInteger);
         break;
       case Unsigned_integer_8_bit:
         return(new Bit8UnsignedInteger);
         break;
       case Signed_integer_16_bit:
         return(new Bit16SignedInteger);
         break;
       case Unsigned_integer_16_bit:
         return(new Bit16UnsignedInteger);
         break;
       case Signed_integer_32_bit:
         return(new Bit32SignedInteger);
         break;
       case Unsigned_integer_32_bit:
         return(new Bit32UnsignedInteger);
         break;
       default:
         return(NULL);
    }
  return(NULL);
}

void
ArrayElem::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
ArrayElem::Print() const
{
}

void
ArrayElem::Delete() const
{
}

void
CellArray::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
CellArray::Print() const
{
  //cout << "CellArray: "<< variable_name << "\n";
}

void
CellArray::Delete() const
{
  //cout << "CellArray: "<< variable_name << "\n";
}

void
CellArray::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  SimpleElem* simple;
  ReadArrayDimension(InBuff, pBuff);
  Cell_number = 1;
  for(unsigned int i=0;i<dimension.size();i++)
    Cell_number *= dimension[i];
  ReadArrayName(InBuff, pBuff);
  flag.character=true;
  for(int i=0;i<Cell_number;i++)
    {
      simple=new SimpleElem;
      Data_Header_t data_header=simple->ReadDataHeader(InBuff, pBuff);
      VCell.push_back(simple);
      simple->DataProceed(data_header, InBuff, pBuff, flag);
    }
}


void
Structure::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
  if(name==variable_name || found)
    {
      found = true;
      vector<string>::const_iterator it2=Structure_Elem_name.begin();
      for(vector<PSimpleElem>::const_iterator it=VCell.begin(); it!=VCell.end(); it++)
        {
          collect_struct.tmp_name = *it2;
          it2++;
          (*it)->Collect(name, found, collect_struct);
        }
    }
}

void
Structure::Print() const
{
  cout << "Structure: " << variable_name << "\n";
  vector<string>::const_iterator it2=Structure_Elem_name.begin();
  int i=0;
  for(vector<PSimpleElem>::const_iterator it=VCell.begin(); it!=VCell.end(); it++)
    {
      cout << ++i << " -> " << *it2 << " :";
      it2++;
      (*it)->Print();
      cout << "\n";
    }
}

void
Structure::Delete() const
{
  vector<string>::const_iterator it2=Structure_Elem_name.begin();
  for(vector<PSimpleElem>::const_iterator it=VCell.begin(); it!=VCell.end(); it++)
    {
      (*it)->Delete();
      delete *it;
    }
}

void
Structure::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  SimpleElem* simple;
  int i;
  ReadArrayDimension(InBuff, pBuff);
  ReadArrayName(InBuff, pBuff);
  ReadStructureNames(InBuff, pBuff);
  flag.no_name=true;
  for(i=0;i<Structure_number;i++)
    {
      simple=new SimpleElem;
      Data_Header_t data_header=simple->ReadDataHeader(InBuff, pBuff);
      simple->DataProceed(data_header, InBuff, pBuff, flag);
      data_header=simple->ReadDataHeader(InBuff, pBuff);
      VCell.push_back(simple);
      simple->DataProceed(data_header, InBuff, pBuff, flag);
    }
}

void
Object::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}


void
Object::Print() const
{
}

void
Object::Delete() const
{
}


void
Object::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  cerr << "Error: Object not implemented\n";
  exit(EXIT_FAILURE);
}

void
CharacterArray::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
  if(name==variable_name || found)
    {
      found = true;
      collect_struct.tmp_name=variable_name;
      vector<PSimpleElem>::const_iterator it=VCell.begin();
      (*it)->Collect(name, found, collect_struct);
    }
}

void
CharacterArray::Print() const
{
  cout << "CharacterArray: " << variable_name << "\n";
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Print();
  cout << "\n";
}

void
CharacterArray::Delete() const
{
  //cout << "CharacterArray: " << variable_name << "\n";
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Delete();
  delete *it;
  //cout << "\n";
}


void
CharacterArray::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  ReadArrayDimension(InBuff, pBuff);
  Matrix_Elem_number = 1;
  for(unsigned int i=0;i<dimension.size();i++)
    {
      Matrix_Elem_number *= dimension[i];
    }
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
SparseArray::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
SparseArray::Print() const
{
  cout << "Sparse Array: " << variable_name << "\n";
}

void
SparseArray::Delete() const
{
  //cout << "Sparse Array: " << variable_name << "\n";
}


void
SparseArray::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  cerr << "Error: Sparse Array not implemented\n";
  exit(EXIT_FAILURE);
}

void
DoublePrecisionArray::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
  if(name==variable_name || found)
    {
      found = true;
      collect_struct.tmp_name=variable_name;
      vector<PSimpleElem>::const_iterator it=VCell.begin();
      (*it)->Collect(name, found, collect_struct);
    }
}


void
DoublePrecisionArray::Print() const
{
  cout << "DoublePrecisionArray: " << variable_name << "\n";
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Print();
  cout << "\n";
}

void
DoublePrecisionArray::Delete() const
{
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Delete();
  delete *it;
}

void
DoublePrecisionArray::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  ReadArrayDimension(InBuff, pBuff);
  Matrix_Elem_number = 1;
  for(unsigned int i=0;i<dimension.size();i++)
    {
      Matrix_Elem_number *= dimension[i];
    }
  if(!flag.no_name)
    ReadArrayName(InBuff, pBuff);
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
SinglePrecisionArray::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
  if(name==variable_name || found)
    {
      found = true;
      collect_struct.tmp_name=variable_name;
      vector<PSimpleElem>::const_iterator it=VCell.begin();
      (*it)->Collect(name, found, collect_struct);
    }
}


void
SinglePrecisionArray::Print() const
{
  cout << "SinglePrecisionArray: " << variable_name << "\n";
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Print();
  cout << "\n";
}

void
SinglePrecisionArray::Delete() const
{
  vector<PSimpleElem>::const_iterator it=VCell.begin();
  (*it)->Delete();
  delete *it;
}

void
SinglePrecisionArray::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  ReadArrayDimension(InBuff, pBuff);
  Matrix_Elem_number = 1;
  for(unsigned int i=0;i<dimension.size();i++)
    {
      Matrix_Elem_number *= dimension[i];
    }
  if(!flag.no_name)
    ReadArrayName(InBuff, pBuff);
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit8SignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}



void
Bit8SignedInteger::Print() const
{
  //cout << "Bit8SignedInteger: \n";
}

void
Bit8SignedInteger::Delete() const
{
  //cout << "Bit8SignedInteger: \n";
}

void
Bit8SignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit8UnsignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}


void
Bit8UnsignedInteger::Print() const
{
  //cout << "Bit8UnsignedInteger: \n";
}

void
Bit8UnsignedInteger::Delete() const
{
  //cout << "Bit8UnsignedInteger: \n";
}

void
Bit8UnsignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit16SignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
Bit16SignedInteger::Print() const
{
  //cout << "Bit16SignedInteger: \n";
}

void
Bit16SignedInteger::Delete() const
{
  //cout << "Bit16SignedInteger: \n";
}

void
Bit16SignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit16UnsignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
Bit16UnsignedInteger::Print() const
{
  //cout << "Bit16UnsignedInteger: \n";
}

void
Bit16UnsignedInteger::Delete() const
{
  //cout << "Bit16UnsignedInteger: \n";
}

void
Bit16UnsignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit32SignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
Bit32SignedInteger::Print() const
{
  //cout << "Bit32SignedInteger: \n";
}

void
Bit32SignedInteger::Delete() const
{
  //cout << "Bit32SignedInteger: \n";
}

void
Bit32SignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

void
Bit32UnsignedInteger::Collect(const string &name, bool found, CollectStruct &collect_struct) const
{
}

void
Bit32UnsignedInteger::Print() const
{
  //cout << "Bit32UnsignedInteger: \n";
}

void
Bit32UnsignedInteger::Delete() const
{
  //cout << "Bit32UnsignedInteger: \n";
}

void
Bit32UnsignedInteger::ReadArray(char* InBuff, int* pBuff, FlagStructure_t flag)
{
  Data_Header_t data_header;
  SimpleElem* simple;
  simple=new SimpleElem;
  VCell.push_back(simple);
  data_header=simple->ReadDataHeader(InBuff, pBuff);
  simple->DataProceed(data_header, InBuff, pBuff, flag);
}

Data_Header_t
MatlabFile::ReadDataHeader(ifstream &MatFile)
{
  Data_Header_t data_header;
  MatFile.read(reinterpret_cast<char*>(&data_header.DataType),sizeof(data_header.DataType));
  MatFile.read(reinterpret_cast<char*>(&data_header.S_Number_of_Bytes),sizeof(data_header.S_Number_of_Bytes));
  if(data_header.S_Number_of_Bytes!=0)
    data_header.Number_of_Bytes=data_header.S_Number_of_Bytes;
  else
    {
      MatFile.read(reinterpret_cast<char*>(&data_header.Number_of_Bytes),sizeof(data_header.Number_of_Bytes));
    }
  if(data_header.Number_of_Bytes<8)
    data_header.Number_of_Bytes = 8;
  return(data_header);
}





void
MatlabFile::MatFileRead(string filename)
{
  ifstream MatFile;
  Data_Header_t data_header;
  SimpleElem *simpl;
  int pBuff;
  FlagStructure_t flag;
  //ArrayElem elem;
  MatFile.open(filename.c_str(),std::ios::in | std::ios::binary);
  if (!MatFile.is_open())
    {
      cerr << filename.c_str() << " Cannot be opened\n";
      exit(EXIT_FAILURE);
    }
  // Read the Header of the Mat-File
  MatFile.read(reinterpret_cast<char*>(&header),sizeof(header));
  do
    {
      data_header=ReadDataHeader(MatFile);
      char* InBuff;
      InBuff = (char*)malloc(data_header.Number_of_Bytes+1);
      MatFile.read(InBuff,data_header.Number_of_Bytes+1);
      pBuff = 0;
      simpl = new SimpleElem;
      VSimpl.push_back(simpl);
      flag.no_name=false;
      flag.character=false;
      simpl->DataProceed(data_header, InBuff, &pBuff, flag);
      free(InBuff);
    }
  while(!MatFile.eof());
  MatFile.close();
}

bool
MatlabFile::Collect(const string &name, CollectStruct &collect_struct) const
{
  for(vector<PSimpleElem>::const_iterator it=VSimpl.begin(); it!=VSimpl.end(); it++)
      (*it)->Collect(name, false, collect_struct);
  return(!(collect_struct.variable_double_name.empty() and collect_struct.variable_string_name.empty()));
}


void
MatlabFile::MatFilePrint()
{
  for(vector<PSimpleElem>::iterator it=VSimpl.begin(); it!=VSimpl.end(); it++)
      (*it)->Print();
}


void
MatlabFile::Delete()
{
  for(vector<PSimpleElem>::iterator it=VSimpl.begin(); it!=VSimpl.end(); it++)
    {
      (*it)->Delete();
      delete *it;
    }
}

/*
int
main(int argc, char** argv)
{
  CollectStruct collect_struct;
  MatlabFile matlab_file;
  matlab_file.MatFileRead("gimf_steady.mat");
  //matlab_file.MatFileRead("essai.mat");
  matlab_file.MatFilePrint();
  bool tmp_b=false;
  //string tmp_s("stored_values");
  tmp_b=matlab_file.Collect("stored_values", collect_struct);
  if(tmp_b)
    {
      int i=0;
      for(map<string,vector<double> >::iterator it2=collect_struct.variable_double_name.begin();it2!=collect_struct.variable_double_name.end();it2++)
        {
          i++;
          cout << i << " " << it2->first.c_str() << " : ";
          for(vector<double>::iterator it=it2->second.begin();it!=it2->second.end();it++)
            cout << *it;
          cout << "\n";
        }
    }
}


*/
