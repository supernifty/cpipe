/*****************************************************************************
  NewGenomeFile.cpp

  (c) 2009 - Aaron Quinlan
  Hall Laboratory
  Department of Biochemistry and Molecular Genetics
  University of Virginia
  aaronquinlan@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/
#include "NewGenomeFile.h"
#include "ParseTools.h"

NewGenomeFile::NewGenomeFile(const QuickString &genomeFilename)
: _maxId(-1)
{
    _genomeFileName = genomeFilename;
    loadGenomeFileIntoMap();
}

NewGenomeFile::NewGenomeFile(const BamTools::RefVector &refVector)
: _maxId(-1)
{
    for (size_t i = 0; i < refVector.size(); ++i) {
        QuickString chrom = refVector[i].RefName;
        CHRPOS length = refVector[i].RefLength;
        _maxId++;
        _chromSizeIds[chrom] = pair<CHRPOS, CHRPOS>(length, _maxId);
    }
}

// Destructor
NewGenomeFile::~NewGenomeFile(void) {
}


void NewGenomeFile::loadGenomeFileIntoMap() {


	ifstream genFile(_genomeFileName.c_str());
	if (!genFile.good()) {
		cerr << "Error: Can't open genome file" << _genomeFileName << "Exiting..." << endl;
		exit(1);
	}
	string sLine;
	vector<QuickString> fields;
	CHRPOS chrSize = 0;
	QuickString chrName;
	while (!genFile.eof()) {
		sLine.clear();
		fields.clear();
		chrSize = 0;
		chrName.clear();
		getline(genFile, sLine);
		Tokenize(sLine.c_str(), fields);
		if (fields.size() != 2) {
			continue;
		}
		chrName = fields[0];
		chrSize = str2chrPos(fields[1]);
		_maxId++;
		_chromSizeIds[chrName] = pair<CHRPOS, int>(chrSize, _maxId);
		_startOffsets.push_back(_genomeLength);
		_genomeLength += chrSize;
		_chromList.push_back(chrName);
	}
	if (_maxId == -1) {
		cerr << "Error: The genome file " << _genomeFileName << " has no valid entries. Exiting." << endl;
		exit(1);
	}
	_startOffsets.push_back(_genomeLength); //insert the final length as the last element
	//to help with the lower_bound call in the projectOnGenome method.
	genFile.close();
}

bool NewGenomeFile::projectOnGenome(CHRPOS genome_pos, QuickString &chrom, CHRPOS &start) {
    // search for the chrom that the position belongs on.
    // add 1 to genome position b/c of zero-based, half open.
    vector<CHRPOS>::const_iterator low =
        lower_bound(_startOffsets.begin(), _startOffsets.end(), genome_pos + 1);
    
    // use the iterator to identify the appropriate index 
    // into the chrom name and start vectors
    CHRPOS i = CHRPOS(low-_startOffsets.begin());
    if (i < 0 || i >= _chromList.size()) {
    	return false; //position not on genome
    }
    chrom = _chromList[i - 1];
    start = genome_pos - _startOffsets[i - 1];
    return true;
}
    
CHRPOS NewGenomeFile::getChromSize(const QuickString &chrom) {
	if (chrom == _currChromName) {
		return _currChromSize;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromSize;
    }
    cerr << "Error: chromosome " << chrom << " is not in the genome file " << _genomeFileName << ". Exiting." << endl;
    return INT_MAX;
}

CHRPOS NewGenomeFile::getChromId(const QuickString &chrom) {
	if (chrom == _currChromName) {
		return _currChromId;
	}
	lookupType::const_iterator iter= _chromSizeIds.find(chrom);
    if (iter != _chromSizeIds.end()) {
    	_currChromName = iter->first;
    	_currChromSize = iter->second.first;
    	_currChromId = iter->second.second;
    	return _currChromId;
    } else {
    	cerr << "Error: requested chromosome " << chrom << " does not exist in the genome file " << _genomeFileName << ". Exiting." << endl;
    	exit(1);
    }
    return INT_MAX;
}
