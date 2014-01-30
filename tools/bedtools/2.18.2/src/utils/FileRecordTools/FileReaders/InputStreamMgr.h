/*
 * InputStreamMgr.h
 *
 *  Created on: Mar 21, 2013
 *      Author: nek3d
 */

#ifndef INPUTSTREAMMGR_H_
#define INPUTSTREAMMGR_H_

using namespace std;

#include "PushBackStreamBuf.h"
#include "InflateStreamBuf.h"
#include "QuickString.h"
#include "api/BamReader.h"
#include "api/internal/io/BgzfStream_p.h"

#include <iostream>

class InputStreamMgr {
public:
	 //Contsructor: 1st arg can be "-" for stdin. set 2nd arg false if fileType already known.
	InputStreamMgr(const QuickString &filename, bool buildScanBuffer = true);
	~InputStreamMgr();
	bool init();
	int read(char *data, size_t dataSize);

	//use getScanBuffer for auto-detection of file types.
//	istream *getFinalStream() { return _finalInputStream; }
	const BTlist<int> &getScanBuffer() const { return _scanBuffer; }
	int getBufferLength() const { return _numBytesInBuffer; }
	void populateScanBuffer();
	const QuickString &getSavedData() const { return _saveDataStr; }
	bool isGzipped() const { return _isGzipped; }
	PushBackStreamBuf *getPushBackStreamBuf() const {return _pushBackStreamBuf; }
//	void getSavedData(QuickString &str) const { str = _saveDataStr; }
	bool isBam() const { return _isBam; }
	BamTools::BamReader *getBamReader() { return _bamReader; }
	bool resetStream();


private:
	QuickString _filename;
	PushBackStreamBuf *_pushBackStreamBuf;
	ifstream *_inputFileStream;
	BTlist<int> _scanBuffer;
	QuickString _saveDataStr;
	InflateStreamBuf *_infStreamBuf;
	istream * _finalInputStream;
	istream *_oldInputStream;
	bool _isStdin;
	bool _isGzipped;
	bool _isBam;
	bool _isBgzipped;
	char *_tmpZipBuf;
	bool _bamRuledOut;
	bool _streamFinished;
	vector<int> _possibleBamCode;
	static const int SCAN_BUFFER_SIZE = 4096; // 4 K buffer
	static const int BAM_SCAN_BUFFER_SIZE = 32768; // 32K
	static const int MIN_SCAN_BUFFER_SIZE = 2048;
	int _numBytesInBuffer; //this will hold the length of the buffer after the scan.
	BamTools::BamReader *_bamReader;
	BamTools::Internal::BgzfStream *_bgStream;

	static const char *FIFO_STRING_LITERAL;
	void readZipChunk();
	bool detectBamOrBgzip(int &numChars, int currChar);
//	void decompressBuffer();

};

#endif /* INPUTSTREAMMGR_H_ */
