// -*- mode: c++ -*-
// $Id: XMLStore.cc 2075 2011-03-03 23:40:05Z axel $

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE /* glibc2 needs this for strptime */
#endif
#ifndef _BSD_SOURCE
#define _BSD_SOURCE /* so we can set TZ in struct tm */
#endif

#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/sax/InputSource.hpp>
#include <xercesc/util/BinInputStream.hpp>
#include <stdlib.h> // for fopen/fclose
#include <limits.h>
#include <bzlib.h>
#include <sys/stat.h>

#include "XMLStore.h"
#include "error.h"

XERCES_CPP_NAMESPACE_USE

// reading from pipe:
class pipeInputStream : public BinInputStream {
  FILE * the_pipe;
public:
  pipeInputStream(const std::string & command):the_pipe(popen(command.c_str(),"r")){
    if(the_pipe==0) FATAL_ERROR("could not open the pipe '" << command << "'");
  }
#if XERCES_VERSION_MAJOR < 3
  virtual unsigned int curPos() const
#else
  virtual XMLFilePos curPos() const
#endif
  {
    return ftell(the_pipe)/sizeof(XMLByte);
  }
#if XERCES_VERSION_MAJOR < 3
  virtual unsigned readBytes(XMLByte* const toFill, const unsigned maxToRead)
#else
  virtual XMLSize_t readBytes(XMLByte* const toFill, const XMLSize_t maxToRead)
#endif
  {
    return fread(toFill, sizeof(XMLByte), maxToRead, the_pipe);
  };
  const XMLCh* getContentType() const{return 0;};
  virtual ~pipeInputStream(){
    fclose(the_pipe);
  };
};

class pipeInputSource : public InputSource {
  std::string the_pipe_name;
public:
  pipeInputSource(const char * pipecommand):the_pipe_name(pipecommand){};
  virtual ~pipeInputSource(){
  };
  virtual BinInputStream * makeStream() const{
    BinInputStream * stream= new pipeInputStream(the_pipe_name);
    if(stream==0) FATAL_ERROR("allocating memory");
    return stream;
  };
};

// writing to pipe:
class pipeFormatTarget : public XMLFormatTarget {
  FILE * the_pipe;
public:
  pipeFormatTarget(const char * pipecommand):
    the_pipe(popen(pipecommand,"w")){
    if(the_pipe==0) 
      FATAL_ERROR("could not open the pipe '" << pipecommand << "'");
  }
  virtual ~pipeFormatTarget(){
    fclose(the_pipe);
  };
#if XERCES_VERSION_MAJOR < 3
  virtual void writeChars(const XMLByte* const toWrite
			  , const unsigned count
			  , XMLFormatter * const )
#else
  virtual void writeChars(const XMLByte* const toWrite
			  , const XMLSize_t count
			  , XMLFormatter * const )
#endif
  {
    if(count!=fwrite(toWrite, sizeof(XMLByte),  count,  the_pipe))
      FATAL_ERROR("trouble writing to pipe");
    return;
  };
  virtual void flushBuffer(){
    fflush(the_pipe);
  };
};

// reading from a bz2 compressed file:
class bz2InputStream : public BinInputStream {
  FILE * f;
  BZFILE * bzf;
  unsigned long int char_position;
  bool eof;
  std::string fname;
public:
  bz2InputStream(const std::string & filename):
    char_position(0),
    eof(false),
    fname(filename)
  {
    // First, set mode to unreadable:
    f=fopen(filename.c_str(),"r");
    if(f==0) FATAL_ERROR("could not open the file '" << fname << "'");
    int bzerror, verbosity=0, small=0;
    bzf=BZ2_bzReadOpen(&bzerror,f,verbosity,small,NULL,0);
    if(bzerror!=BZ_OK){
      FATAL_ERROR("trouble preparing '" << fname << "' for compression");
    }
//     REPORT(bzf);
//     WARNING("opened");
  }
#if XERCES_VERSION_MAJOR < 3
  virtual unsigned int curPos() const
#else
  virtual XMLFilePos curPos() const
#endif
  {
    return char_position/sizeof(XMLByte);
  }
#if XERCES_VERSION_MAJOR < 3
  virtual unsigned readBytes(XMLByte* const toFill, const unsigned maxToRead)
#else
  virtual XMLSize_t readBytes(XMLByte *const toFill,
				 const XMLSize_t maxToRead)
#endif
    {
//     REPORT(eof);
    if(eof) return 0;
    int bzerror;
//     REPORT(bzf);
    int chars_read=
      BZ2_bzRead(&bzerror, bzf, toFill, sizeof(XMLByte)*maxToRead);
    ERROR_TEST(bzerror,BZ_PARAM_ERROR);
    ERROR_TEST(bzerror,BZ_SEQUENCE_ERROR);
    ERROR_TEST(bzerror,BZ_IO_ERROR);
    ERROR_TEST(bzerror,BZ_UNEXPECTED_EOF);
    ERROR_TEST(bzerror,BZ_DATA_ERROR);
    ERROR_TEST(bzerror,BZ_DATA_ERROR_MAGIC);
    ERROR_TEST(bzerror,BZ_MEM_ERROR);
    if(bzerror!=BZ_OK && bzerror!=BZ_STREAM_END)
      FATAL_ERROR("trouble reading '" << fname << "'");
    if(bzerror==BZ_STREAM_END)
      eof=true;
    if( (chars_read%sizeof(XMLByte)) != 0 )
      FATAL_ERROR("could not read a full XMLByte");
    char_position+=chars_read;
    return chars_read/sizeof(XMLByte);
  };
  const XMLCh* getContentType() const{return 0;};
  virtual ~bz2InputStream(){
    int bzerror;
    BZ2_bzReadClose(&bzerror, bzf);
    if(bzerror!=BZ_OK)
      FATAL_ERROR("trouble closing '" << fname << "'");
    if(fclose(f)){
      perror("");
      FATAL_ERROR("trouble closing '" << fname << "'");
    }
//     REPORT(bzf);
//     WARNING("closed");
//     WARNING("bz2InputStream destroyed OK");
  };
};

class bz2InputSource : public InputSource {
  std::string the_filename;
  mutable bool produced_stream;
public:
  bz2InputSource(const char * filename):
    the_filename(filename),
    produced_stream(false)
  {};
  virtual ~bz2InputSource(){
    if(produced_stream)
      ;//WARNING("bz2InputSource destroyed OK");
    else
      WARNING("bz2InputSource destroyed before producing stream");
  };
  virtual BinInputStream * makeStream() const{
    BinInputStream * stream= new bz2InputStream(the_filename);
    if(stream==0) FATAL_ERROR("allocating memory");
    produced_stream=true;
    return stream;
  };
};

// writing with bzip2 compression:
class bz2FormatTarget : public XMLFormatTarget {
  FILE * f;
  BZFILE * bzf;
  std::string fname;
public:
  bz2FormatTarget(const char * filename):fname(filename){
    // First, set mode to unreadable:
    mode_t old_umask=umask(S_IRUSR | S_IRGRP | S_IROTH | S_IWGRP | S_IWOTH);
    f=fopen(filename,"w");
    if(f==0) FATAL_ERROR("could not open the file '" << fname << "'");
    umask(old_umask); // reset old umask

    int bzerror,verbosity=0,blockSize100k=9,workFactor=0/*=use default*/;
    bzf=BZ2_bzWriteOpen(&bzerror,f,blockSize100k,
			verbosity,workFactor);
    if(bzf==0)
      FATAL_ERROR("could not initialize compresssion of '" 
		  << fname << "'.");
  }
  virtual ~bz2FormatTarget(){
    int bzerror;
    BZ2_bzWriteClose(&bzerror,bzf,0,NULL,NULL);
    ERROR_TEST(bzerror,BZ_SEQUENCE_ERROR);
    ERROR_TEST(bzerror,BZ_IO_ERROR);
    if(bzerror!=BZ_OK){
      FATAL_ERROR("trouble closing '" << fname << "'");
    }
    if(fclose(f)){
      perror("");
      FATAL_ERROR("trouble closing '" << fname << "'");
    }
    chmod(fname.c_str(),S_IRUSR|S_IWUSR|S_IRGRP);

  };
#if XERCES_VERSION_MAJOR < 3
  virtual void writeChars(const XMLByte* const toWrite
			  , const unsigned count
			  , XMLFormatter * const formatter)
#else
  virtual void writeChars(const XMLByte* const toWrite
			  , const XMLSize_t count
			  , XMLFormatter * const formatter)
#endif
  {
    if((long int)(sizeof(XMLByte)/sizeof(char))*count>INT_MAX){
      unsigned int half=count/2;
      writeChars(toWrite,half,formatter);
      writeChars(toWrite+half,count-half,formatter);
    }else{
      int bzerror;
      BZ2_bzWrite(&bzerror,bzf,const_cast<XMLByte *>(toWrite),
		  (sizeof(XMLByte)/sizeof(char))*count);
      if(bzerror != BZ_OK)
	FATAL_ERROR("trouble writing to '" << fname << "'");
    }
    return;
  };
  virtual void flushBuffer(){
    // Flushing makes little sense for compressed files.  We could do
    // a close-open combination, but probably this will not give what
    // we want.
  };
};

const char * XMLStore::the_implementation_descriptor;

XMLStore::XMLStore(){
  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    FATAL_ERROR("could not initialize XERCES");
  }
  the_implementation = 
    DOMImplementationRegistry::
    getDOMImplementation(toXMLCh(the_implementation_descriptor));
  if(!the_implementation) throw XMLStoreException("trouble registring at XERCES");
  the_doc = the_implementation->createDocument(0, toXMLCh("FoodWebData"), 0);
  if(!the_doc) throw XMLStoreException("trouble registring at XERCES");
  parser=0;
}

XMLStore::XMLStore(const std::string & fn){
  const char * filename=fn.c_str();

  // check if file readable
  FILE * fd;
  if(!(fd=fopen(filename,"r"))){
    FATAL_ERROR("cannot read file \"" << filename << '"');
  }else{
    fclose(fd);
  }

  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    FATAL_ERROR("could not initialize XERCES");
  }
  parser = new XercesDOMParser;  //parser is private in XMLStore
  parser->useImplementation(toXMLCh(the_implementation_descriptor));
  parser->setValidationScheme(XercesDOMParser::Val_Auto);
  parser->setDoNamespaces(false);
  parser->setDoSchema(false);
  parser->setValidationSchemaFullChecking(false);
  //parser->setCreateEntityReferenceNodes(); // leave at default
  
  try{
    if(strlen(filename)>3 && 0==strcmp(filename+strlen(filename)-3,".gz")){
      const char * gzip = "zcat ";
      char * pipecommand = new char[strlen(filename)+strlen(gzip)+1];
      pipecommand[0]=0;
      strcat(pipecommand,gzip);
      strcat(pipecommand,filename);
      pipeInputSource pis(pipecommand);
      parser->parse(pis);
      delete[] pipecommand;
    }else if(strlen(filename)>4 && 0==strcmp(filename+strlen(filename)-4,".bz2")){
#if 0 // dead code
      const char * gzip = "bzcat ";
      char * pipecommand = new char[strlen(filename)+strlen(gzip)+1];
      pipecommand[0]=0;
      strcat(pipecommand,gzip);
      strcat(pipecommand,filename);
      pipeInputSource pis(pipecommand);
      parser->parse(pis);
      delete[] pipecommand;
#else
      parser->parse(bz2InputSource(filename));
#endif
    }else{
      parser->parse(filename);
    }
  }
  catch (const OutOfMemoryException&)
    {
      FATAL_ERROR("out of memmory");
    }
  catch (const XMLException& e)
    {
      FATAL_ERROR(tochar(e.getMessage()) << "(during parsing)");
    }
  catch (const DOMException& e)
    {
      const unsigned int maxChars = 2047;
      XMLCh errText[maxChars + 1];
      
      if (!DOMImplementation::loadDOMExceptionMsg(e.code, errText, maxChars))
	errText[0]=0;
      FATAL_ERROR(tochar(errText) << " (error code " << e.code << ")");
    }
  catch (...){
    FATAL_ERROR("some other error");
  }
  
  the_doc = parser->getDocument();
  // deleting here would cause a segfault:
  //delete parser;  // not deleting here causes a memory leak.
  the_implementation = 0; //???? how can I get the parsers implementation??
}

XMLStore::~XMLStore(){
  //the_doc->release();  // this and delete parser together causes a problem
  if(parser){
    parser->resetDocumentPool();
    delete parser; 
  }else{
    the_doc->release();
  }
  XMLPlatformUtils::Terminate();
}

XMLStore & XMLStore::put(permanent *data, const char * tag){
  if(remember(the_doc,the_doc->getDocumentElement(),rememberWrite)
     .sync(tag,data)){
    throw XMLStoreException("problem while putting data");
  }
  return *this;
}

XMLStore & XMLStore::get(permanent *data, const char * tag){
  if(remember(the_doc,the_doc->getDocumentElement(),rememberRead)
     .sync(tag,data)){
    WARNING("problem while getting " << tag);
    throw XMLStoreException("problem while getting data");
  }
  return *this;
}

XMLStore & XMLStore::save(const char * filename){
  try
    {
#if XERCES_VERSION_MAJOR < 3
      DOMWriter *theSerializer = the_implementation->createDOMWriter();
      // set user specified output encoding
      theSerializer->setEncoding(0/*gOutputEncoding*/); //default (UTF-8)
      if(theSerializer->canSetFeature(XMLUni::fgDOMWRTSplitCdataSections,true))
	theSerializer->setFeature(XMLUni::fgDOMWRTSplitCdataSections,true);
      if(theSerializer->canSetFeature(XMLUni::fgDOMWRTDiscardDefaultContent,true))
	theSerializer->setFeature(XMLUni::fgDOMWRTDiscardDefaultContent,true);
      if(theSerializer->canSetFeature(XMLUni::fgDOMWRTFormatPrettyPrint,true))
	theSerializer->setFeature(XMLUni::fgDOMWRTFormatPrettyPrint,true);
#else
      DOMLSSerializer *theSerializer = the_implementation->createLSSerializer();
      theSerializer->getDomConfig()
	->setParameter(toXMLCh("format-pretty-print"),true);
#endif
	
      //
      // Plug in a format target to receive the resultant
      // XML stream from the serializer.
      //
      // StdOutFormatTarget prints the resultant XML stream
      // to stdout once it receives any thing from the serializer.
      //
      XMLFormatTarget *myFormTarget;
      if(strlen(filename)>3 && 0==strcmp(filename+strlen(filename)-3,".gz")){
	const char * gzip = "umask 477;gzip --best > ";
	const char * unlock = ";chmod 640 ";
	char * pipecommand = 
	  new char[2*strlen(filename)+strlen(gzip)+strlen(unlock)+1];
	pipecommand[0]=0;
	strcat(pipecommand,gzip);
	strcat(pipecommand,filename);
	strcat(pipecommand,unlock);
	strcat(pipecommand,filename);
	myFormTarget = new pipeFormatTarget(pipecommand);
	delete[] pipecommand;
      }else      
	if(strlen(filename)>4 && 0==strcmp(filename+strlen(filename)-4,".bz2")){
#if 0 // dead code
	const char * gzip = "umask 477;bzip2 --best > ";
	const char * unlock = ";chmod 640 ";
	char * pipecommand = 
	  new char[2*strlen(filename)+strlen(gzip)+strlen(unlock)+1];
	pipecommand[0]=0;
	strcat(pipecommand,gzip);
	strcat(pipecommand,filename);
	strcat(pipecommand,unlock);
	strcat(pipecommand,filename);
	myFormTarget = new pipeFormatTarget(pipecommand);
	delete[] pipecommand;
#else
	myFormTarget = new bz2FormatTarget(filename);
#endif
      }else{
	myFormTarget = new LocalFileFormatTarget(filename);
      }
      //
#if XERCES_VERSION_MAJOR < 3
      theSerializer->writeNode(myFormTarget, *the_doc);
#else
      DOMLSOutput * myOutput = the_implementation -> createLSOutput();
      myOutput->setByteStream(myFormTarget);
      theSerializer->write(the_doc, myOutput);
      delete myOutput;
#endif	
      delete theSerializer;
      //
      // Filter, formatTarget and error handler
      // are NOT owned by the serializer.
      //
      delete myFormTarget;
	
    }
  catch (const OutOfMemoryException&)
    {
      FATAL_ERROR("out of memory");
    }
  catch (XMLException& e)
    {
      FATAL_ERROR(tochar(e.getMessage()) << "(during parsing)");
    }

  return *this;
}  

int XMLStore::node_count(const char * tag) {
  DOMNodeList*    nodeList = the_doc->getElementsByTagName(toXMLCh(tag));
  return nodeList->getLength();
};

bool XMLStore::older_than(const char * bug_correction_date){
  // used to make code compatibile with buggy older versions:
  if(old_compilation_time.size()==0){
    try{
      this->get(old_compilation_time,"compilation_time");
      REPORT(old_compilation_time);
    }catch(XMLStoreException){
      WARNING("could not find compilation_time in XML data");
    } 
  }
    
  //typical date: May 23 2006, 14:14:16
  struct tm time_data;
  char * strptime_retval;
  strptime_retval=
    strptime(old_compilation_time.c_str(),"%b %d %Y, %T",&time_data);
  if(old_compilation_time.size()>0 && !strptime_retval) 
    FATAL_ERROR("could not parse compliation time of input web");
  //and now the time of the bug fix: 2006/08/01 05:45:21 UTC
  struct tm bug_fix_time;
  strptime_retval=
    strptime(bug_correction_date,"%b %d %Y, %T",&bug_fix_time);
  if(!strptime_retval) 
    FATAL_ERROR("could not parse bug correction date");
//   REPORT(mktime(&time_data));
//   REPORT(mktime(&bug_fix_time));
  if(old_compilation_time.size()==0 || 
     mktime(&time_data)<mktime(&bug_fix_time) ){
    return true;
  }else{
    return false;
  }
}


/*
#include "evaluate.h"
#include "random.h"
#include "Statistics.h"
#include "cfgList.h"
#include "attributes.h"

int main(void){
  NewWeb web;
  int i;
  XMLStore("file:///tmp/web.xml").get(i,"someVariable");
  XMLStore("/tmp/web2.xml.gz").get(&web,"FoodWeb");
  set_random_seed(123);
  web.speciate();
  web.speciate();
  try{
    XMLStore().put(&web,"FoodWeb").put(i,"someVariable").save("/tmp/web2.xml.gz");
  }catch(...){
    FATAL_ERROR("trouble saving");
  }
}
*/
