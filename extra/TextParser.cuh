#ifndef _TEXTPARSER_INC
#define _TEXTPARSER_INC

#include <stdio.h>
/*
A classic, used everywhere to read formatted text files without using Lex
and Yacc. Written in C around 1995
Re-vamped a million times, and included in several programs
Converted to C++ in 1.7.98 to use it in the Rayman project.
*/

class TextParser
{
	char *data;
	unsigned int sl;
	unsigned int size;
public:
	TextParser();
	TextParser(const char*);
	~TextParser();

	bool create(const char *);
	char *getword();
	char *getcommaword();
	int getint();
	double getfloat();

	int countword(char *);
	int countwordfromhere(char *);
	int countchar(char);
	void reset();
	void destroy();
	void goback();
	bool seek(const char *);
	int eof();

	// Funcions afegides.
	int CountObjs();
};



#define TextParserSafeCall( err ) __textParserSafeCall( err, __FILE__, __LINE__ )
inline void __textParserSafeCall(bool err, const char *file, const int line)
{

	if (true != err) {
		fprintf(stderr, "TextParser error in file '%s', line %d\n", __FILE__, __LINE__);
		exit(EXIT_FAILURE);
	}


	return;
}





#endif




