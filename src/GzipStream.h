

/*
The MIT License (MIT)

Copyright (c) 2019 Pierre Lindenbaum @yokofakun Institut du Thorax France

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/
#ifndef GZIP_STREAM
#define GZIP_STREAM

#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include <zlib.h>
#include <streambuf>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>
#include <cerrno>

using namespace std;

class GzipStreamBuf : public std::basic_streambuf<char>
	{
	private:
		unsigned int buff_size;
		char* buffer;
		gzFile gzin;

		void checkError()
			{
			int ret = 0;
			const char* msg = ::gzerror(this->gzin,&ret);
			if(ret==0) return;
			if(msg==NULL)
				{
				throw std::runtime_error("GZLIB: I/O error");
				}
			else
				{
				throw std::runtime_error(msg);
				}
			}
	



		
		void _init(const char* fname,unsigned int buff_size)
			{
			assert(buff_size>0);
			this->gzin = ::gzopen(fname,"r");
			if(this->gzin == NULL)
				{
				std::ostringstream msg;
				msg << "File Parser: could not open file: " << fname << std::endl;
				throw std::runtime_error(msg.str());
				}

			this->buff_size=buff_size;
			this->buffer=new char[buff_size];
		
			setg(	(char*)&this->buffer[0],
				(char*)&this->buffer[this->buff_size],
				(char*)&this->buffer[this->buff_size]
				);
			
			}
		
	
	
	public:
		GzipStreamBuf(const char* fname)
			{
			_init(fname,BUFSIZ);
			}
			
		virtual ~GzipStreamBuf()
			{
			if(gzin!=NULL) ::gzclose(this->gzin);
			if(this->buffer!=NULL) delete [] this->buffer;
			}
		
		virtual int underflow ( )
			{
			int nRead =0;
			if(gzeof(this->gzin)) return EOF;
			
			if( ( nRead = ::gzread(this->gzin,this->buffer,this->buff_size) ) <= 0 ) {
				checkError();
				return EOF;
				}

			setg(	(char*)this->buffer,
				(char*)&this->buffer[1],
				(char*)&this->buffer[nRead]
				);
			
			return this->buffer[0];
			}
	};



#endif
