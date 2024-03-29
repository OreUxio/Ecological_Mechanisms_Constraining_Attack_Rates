// -*- mode: c++ -*-
// $Id: remember.h 1431 2009-05-04 12:22:40Z axel $
#ifndef _REMEMBER_H_
#define _REMEMBER_H_

enum remember_t {rememberWrite,rememberRead};

#include <xercesc/dom/DOM.hpp>
#include <sstream>
#include <iomanip>
#include "error.h"

class tochar
{
private :
  char*   fLocalForm;
public :
  tochar(const XMLCh* const toTranscode)
  {
    // Call the private transcoding method
    fLocalForm = XERCES_CPP_NAMESPACE_QUALIFIER XMLString::transcode(toTranscode);
  }
  ~tochar()
  {
    XERCES_CPP_NAMESPACE_QUALIFIER XMLString::release(&fLocalForm);
  }
  operator const char * () const
  {
    return fLocalForm;
  }
};

class toXMLCh
{
private :
  XMLCh*   fLocalForm;
public :
  toXMLCh(const char* const toTranscode)
  {
    // Call the private transcoding method
    fLocalForm = XERCES_CPP_NAMESPACE_QUALIFIER XMLString::transcode(toTranscode);
  }
  ~toXMLCh()
  {
    XERCES_CPP_NAMESPACE_QUALIFIER XMLString::release(&fLocalForm);
  }
  operator const XMLCh * () const
  {
    return fLocalForm;
  }
};

class remember;

/// Base class for permanent classes on which XMLstore operates.
/** The mapping from data to XML tags and entries is done in
    data_mapping(). Often, it is sufficient to call PERMANENT(\a a)
    for all atomic data \a a, and PERMANENTP(& \a c) for classes \a c
    that need to be remembered in data_mapping(). */
class permanent{
public:
  remember * the_remember;
  void fix(remember * rem){
    the_remember=rem;
  }
  virtual ~permanent(){};
  virtual void data_mapping()=0;
};

/// A technicality.
class remember{
private:
  XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument * the_doc;
  XERCES_CPP_NAMESPACE_QUALIFIER DOMElement * the_parent_node;
  remember_t the_read_or_write;
  friend class permanent;
public:
  remember(XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument * doc,
	   XERCES_CPP_NAMESPACE_QUALIFIER DOMElement * parent_node,
	   remember_t rw):
    the_doc(doc),
    the_parent_node(parent_node),
    the_read_or_write(rw){};
  template<typename T>
  int sync(const char * node_name,T & x);
  int sync(const char * node_name,permanent * c);
  int sync(const char * node_name,std::string & x);
  remember_t read_or_write(){return the_read_or_write;};
};
  
template<typename T>
int remember::sync(const char * node_name,T & x){
  XERCES_CPP_NAMESPACE_QUALIFIER DOMNode* child=0;
  switch(read_or_write()){
  case rememberWrite:
    {
      child = the_doc->createElement(toXMLCh(node_name));
      the_parent_node->appendChild(child);
      std::ostringstream os;
      os << std::setprecision(20) << x;
      XERCES_CPP_NAMESPACE_QUALIFIER DOMText * contents;
      contents = the_doc->createTextNode(toXMLCh(os.str().c_str()));
      child->appendChild(contents);
    }
    break;
  case rememberRead:
    {
        //find the relevant data
      XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList * children =
	the_parent_node->getChildNodes();
      for(unsigned int i=0;i<children->getLength();i++){
	//std::cout << tochar(children->item(i)->getNodeName()) << std::endl;
	if(children->item(i)->getNodeType()==XERCES_CPP_NAMESPACE_QUALIFIER DOMNode::ELEMENT_NODE)
	  child=
	    static_cast<XERCES_CPP_NAMESPACE_QUALIFIER DOMElement*>(children->item(i));
	else
	  child=0;
	if(child && 0==strcmp(tochar(child->getNodeName()),node_name))
	  break;
	else
	  child=0;
      }
      if(!child){
	return 1;
      }

      XERCES_CPP_NAMESPACE_QUALIFIER DOMNode * contents=
	child->getFirstChild();
      if(!contents){
	return 1;
      }
      if(contents->getNodeType()!=
	 XERCES_CPP_NAMESPACE_QUALIFIER DOMNode::TEXT_NODE)
	FATAL_ERROR("syntax_error");
      XERCES_CPP_NAMESPACE_QUALIFIER DOMText * text=
	(XERCES_CPP_NAMESPACE_QUALIFIER DOMText *) contents;
      std::istringstream is;
      is.str(std::string(tochar(text->getData())));
      is >> x;
      the_parent_node->removeChild(child);
      child->release();
    }
    break;
  default:
    FATAL_ERROR("wrong remember direction");
  }
  return 0;
};


#define PERMANENT(VAR) permanent::the_remember->sync(#VAR,(VAR))
#define PERMANENTP(VAR) permanent::the_remember->sync(#VAR,&(VAR))

#endif // _REMEMBER_H_
