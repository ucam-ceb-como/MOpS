#ifndef FL_ERROR_HANDLER_H
#define FL_ERROR_HANDLER_H

#include <string>
using namespace std;

namespace FlameLab{
	class ErrorHandler{
	public:
		string errorString;
		int errorNum;

		ErrorHandler(string errStr, int num){
			errorString=errStr;
			errorNum = num;
		}
	};
}

#endif
	

