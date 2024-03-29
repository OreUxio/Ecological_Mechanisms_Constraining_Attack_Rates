// -*- mode: c++ -*-
// $Id: XMLStore.h 1782 2010-04-30 13:25:07Z axel $
#ifndef _XMLSTORE_H_
#define _XMLSTORE_H_

#include <string>
// some of these includes might better go into XMLStore.cc
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include "remember.h"

#if XERCES_VERSION_MAJOR < 3 
typedef long int XMLFilePos;
#endif

/// Used by XMLStore.
class XMLStoreException{
  const char * message;
 public:
  XMLStoreException(const char * s):message(s){};
};

/// Reads or writes specializations of \a permanent object to files in XML.
/** Use put(...) or get(...). */
class XMLStore{
private:
  static const char * the_implementation_descriptor;
#if XERCES_VERSION_MAJOR >= 3
  XERCES_CPP_NAMESPACE_QUALIFIER DOMImplementationLS* the_LSimplementation;
#endif
  XERCES_CPP_NAMESPACE_QUALIFIER DOMImplementation* the_implementation;
  XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument*  the_doc;
  XERCES_CPP_NAMESPACE_QUALIFIER XercesDOMParser* parser;
  std::string old_compilation_time;
public:
  XMLStore(); ///creates a new DOM tree
  XMLStore(const std::string & filename); ///creates a DOM tree, loads XML data
  ~XMLStore(); ///destroyes the DOM tree, without saving (!!)
  XMLStore & save(const char * filename); ///saves XML data
  XMLStore & save(const std::string & filename){
    return save(filename.c_str());
  }; 
  template<typename T> XMLStore & put(T & data, const char * tag);
  template<typename T> XMLStore & get(T & data, const char * tag);
  XMLStore & put(permanent * data, const char * tag);
  XMLStore & get(permanent * data, const char * tag);
  int node_count(const char * tag="*");
  bool older_than(const char * bug_correction_date);
};

template<typename T> XMLStore & XMLStore::put(T & data, const char * tag){
  if(remember(the_doc,the_doc->getDocumentElement(),rememberWrite)
     .sync(tag,data)){
    throw XMLStoreException("problem while putting data");
  }
  return *this;
}

template<typename T> XMLStore & XMLStore::get(T & data, const char * tag){
  if(remember(the_doc,the_doc->getDocumentElement(),rememberRead)
     .sync(tag,data)){
    WARNING("problem while getting " << tag);
    throw XMLStoreException("problem while getting data");
  }
  return *this;
}

// template <typename T>
// class on_destroy_put{
//   T & the_loc;
//   std::string the_tag;
//   XMLStore & the_store;
// public:
//   on_destroy_put(T & x,const char * tag,XMLStore & s):
//     the_loc(x),the_tag(tag),the_store(s){};
//   ~on_destroy_put(){
//     the_store.put(the_loc,the_tag.c_str());
//   }
// };

// class on_destroy_putp{
//   permanent & the_loc;
//   std::string the_tag;
//   XMLStore & the_store;
// public:
//   on_destroy_putp(permanent & x,const char * tag,XMLStore & s):
//     the_loc(x),the_tag(tag),the_store(s){};
//   ~on_destroy_putp(){
//     the_store.put(& the_loc,the_tag.c_str());
//   }
// };

// class on_destroy_save_XMLStore : public XMLStore{
//   std::string filename;
// public:
//   on_destroy_save_XMLStore(char * f):filename(f){};
//   ~on_destroy_save_XMLStore(){
//     save(filename.c_str());
//   }
// };
  

#endif // _XMLSTORE_H_
