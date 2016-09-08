/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2000 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Ch. Godin (christophe.godin@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers    : amldevlp@cirad.fr
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */





#ifndef _SYMBOL_SAMIR_DEF_
#define _SYMBOL_SAMIR_DEF_

#include "definitions.h"
#include "tool/rw_slist.h"
#include "tool/rw_cstring.h"
//#include <rw/tvslist.h>
//#include <rw/cstring.h>

class MySymbol
{
 public:

  // constructor

  MySymbol(const char*,
           const char*,
           char);

  MySymbol(const char);
  MySymbol(const MySymbol&);
  MySymbol();

  // destructor

  ~MySymbol();

  // Selectors

  const RWCString& getSymbolFileName() const;
  const RWCString& getSymbolName() const;
  char getSymbol() const {return _symbol;};
  // Operators

  Boolean operator==(const MySymbol&) const;
  const MySymbol& operator=(const MySymbol&);

  Boolean isValid();

 private:

  RWCString* _fileName;   // contains the geometric representation of the symbol (cf dressing file)
  RWCString* _symbolName;
  char _symbol;

  void clearAndDestroy();
  void symbolCopy(const MySymbol&);

};

typedef RWTValSlist<MySymbol> SymbolsTable;


class SymbLib
{
 public:

  SymbLib();
  ~SymbLib();

  const MySymbol lookForASymbol(const char) const;
  Boolean addSymbol(const MySymbol);
  const SymbolsTable& getSymbolsTable() const;
  void clear();

  void report() const;


 private:

  SymbolsTable* _symbList;

};

#endif




