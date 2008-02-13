#include "stringf.h"

namespace ComoString {
    inline bool rtrimable(std::wstring & str) {
	    if (str.length() == 0) {
		    return false;
	    } else if ((str[str.length()-1] == L' ')||(str[str.length()-1] == L'\n')||(str[str.length()-1] == L'\t')||(str[str.length()-1] == L'\r')) {
		    return true;
	    } else {
		    return false;
	    }
    }
    inline bool ltrimable(std::wstring & str) {
	    if (str.length() == 0) {
		    return false;
	    } else if ((str[0] == L' ')||(str[0] == L'\n')||(str[0] == L'\t')||(str[0] == L'\r')) {
		    return true;
	    } else {
		    return false;
	    }
    }
    std::wstring trim(std::wstring str) {

        while (rtrimable(str)) {
            str.erase(str.length()-1,1);
        }
        while (ltrimable(str)) {
            str.erase(0,1);
        }
        return str;
    }

    std::wstring insertSpace(const int depth) {
        std::wstring space = L"";
	    for (int i=0; i < depth; i++) {
		    space.append(L"    ");
	    }
        return space;
    }
}
