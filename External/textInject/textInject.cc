#include <windows.h>
#include "mex.h"

int main(char* windowNameBuf,POINT coordsInTheEditField,char* textToTypeBuf)

{

    int buflen = strlen(textToTypeBuf);

	HWND hwnd = FindWindow(NULL,TEXT(windowNameBuf));

	if (hwnd != 0)
		for(int i = 0; i < buflen; i++)
			SendMessage(ChildWindowFromPoint(hwnd,coordsInTheEditField),WM_CHAR,textToTypeBuf[i],NULL);

	else
		mexPrintf("Could not find window of that name...\n");
	
	return 0;

}


// The gateway routine
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{

	// First argument is buffer to hold window name
	char *windowNameBuf;
	mwSize buflen1;
	const mxArray *string_array_ptr1 = prhs[0];
	buflen1 = mxGetNumberOfElements(string_array_ptr1) + 1;
	windowNameBuf = (char *) mxCalloc(buflen1, sizeof(char));
	mxGetString(string_array_ptr1, windowNameBuf, buflen1);

	// Next 2 arguments are pixel co-ordinates that lie within the text field to which we would like to inject text.
	// Co-ordinates relative to the window. Try a screenshot and look at co-ordinates in MSPaint if having trouble.
	int x = (int) *(mxGetPr(prhs[1]));
	int y = (int) *(mxGetPr(prhs[2]));
	POINT coordsInTheEditField;
	coordsInTheEditField.x = x;
	coordsInTheEditField.y = y;

	// Last argument is the text to type. If you want escape characters you will need to encapsulate this argument in sprintf()
	// when calling from Matlab
	char *textToTypeBuf;
	mwSize buflen2;
	const mxArray *string_array_ptr2 = prhs[3];
	buflen2 = mxGetNumberOfElements(string_array_ptr2) + 1;
	textToTypeBuf = (char *) mxCalloc(buflen2, sizeof(char));
	mxGetString(string_array_ptr2, textToTypeBuf, buflen2);

	// Call the C subroutine.
	main(windowNameBuf,coordsInTheEditField,textToTypeBuf);

}