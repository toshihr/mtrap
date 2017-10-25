#include <iostream>
#include <string>
#include <fstream>
//#include <vector>
#include <cassert>
#include <cstdlib>
#include "utility.h"
#include "scorematrix.h"
#include "config.h"

#ifndef UNIX
#include <Windows.h>
#endif

using namespace std;

VS fileSearchDirList;

int getIndex(const vec_pString& v, const std::string& s)
{
    int i=0;
    for(vec_pString::const_iterator p=v.begin(); p!=v.end(); ++i,++p)
    {
        if(**p == s) return i;
    }
    return i;
}

void strReplace(std::string& str, const std::string& from, const std::string& to)
{
    std::string::size_type pos = 0;
    while(pos = str.find(from, pos), pos != std::string::npos) {
        str.replace(pos, from.length(), to);
        pos += to.length();
    }
}

void adjustDirDelim(std::string& str)
{
#ifdef UNIX
    strReplace(str,"\\","/");
#else
    strReplace(str,"/","\\");
#endif
}

void initFileSearchDir()
{
    fileSearchDirList.clear();

/*
    // read an environment variable
    char* env = getenv("MTRAP_RESOURCE");
    if(env)
    {
        fileSearchDirList.push_back(string(env) + "/");
    }
*/

#ifdef UNIX
    // set predefined search dir.
    fileSearchDirList.push_back(string(DATA_DIR) + "/");
#else
    // for VC code
    HKEY hkResult;
    LONG lret = RegOpenKeyEx(HKEY_LOCAL_MACHINE, TEXT("Software\\Microsoft\\Windows\\CurrentVersion\\App Paths\\mtrap.exe"),0,KEY_READ,&hkResult);
    if(lret == ERROR_SUCCESS)
    {
        TCHAR lpData[256];
        DWORD dwType;
        DWORD dwSize = sizeof(lpData)/sizeof(lpData[0]);
        if(RegQueryValueEx(hkResult, 0, 0, &dwType, (LPBYTE)lpData, &dwSize) == ERROR_SUCCESS)
        {
            string fullPath = string(lpData);
            string path = fullPath.substr(0,fullPath.size()-string("\\mtrap.exe").size());
            fileSearchDirList.push_back(path + "/");
        }else{
            putError("no MTRAP registry key. please install again.", false);
        }
        RegCloseKey(hkResult);
    }else{
        putError("no registry key. please install MTRAP correctly.",false);
    }
#endif

    // adjust
    for(VS::iterator i=fileSearchDirList.begin(); i!=fileSearchDirList.end(); ++i)
    {
        adjustDirDelim(*i);
    }
}

string fixFilePath(const string& str)
{
    // ---- [FIRST PRIORITY] treat the reserved hard coding files ---
    if(str == "CID" || str == "CPAM"
            || str == "CGONNET40"
            || str == "CGONNET80"
            || str == "CGONNET120"
            || str == "CGONNET160"
            || str == "CGONNET250"
            || str == "CGONNET300"
            || str == "CGONNET350" )
    {
        return str;
    } else {
        // ---- [SECOND PRIORITY] just try open ---
        {
            string fullName = str;
            adjustDirDelim(fullName);
            ifstream ifs(fullName.c_str(), ios::in);
            if(ifs)
            {
                ifs.close();
                return fullName;
            }
        }

        // ---- [THIRD PRIORITY] search directories ---
        for(vector<string>::iterator i=fileSearchDirList.begin(); i!=fileSearchDirList.end(); ++i)
        {
            string fullName = *i + str;
            adjustDirDelim(fullName);
            ifstream ifs(fullName.c_str(), ios::in);
            if(ifs)
            {
                ifs.close();
                return fullName;
            }
        }
        putError("Can not find " + str + " in the following dir.", false);
        for(vector<string>::iterator i=fileSearchDirList.begin(); i!=fileSearchDirList.end(); ++i)putError(*i, false);
        putError("You should use -I option to add a search directory.");
        return "";
    }
}

void putError(string& str, bool finish)
{
    cerr << "ERR: " << str << endl;
    if(finish)exit(1);
}

void putError(const string& str, bool finish)
{
    cerr << "ERR: " << str << endl;
    if(finish)exit(1);
}

void putLog(std::string& str)
{
    cerr << str << endl;
}

void putLog(const std::string& str)
{
    cerr << str << endl;
}

void putWarning(std::string& str)
{
    cerr << "WARNING: " << str << endl;
}

void putWarning(const std::string& str)
{
    cerr << "WARNING: " << str << endl;
}

void mySplit(const string& str, const string& delim, vector<string>& parts) {
  size_t start, end = 0;
  while (end < str.size()) {
    start = end;
    while (start < str.size() && (delim.find(str[start]) != string::npos)) {
      start++;  // skip initial whitespace
    }
    end = start;
    while (end < str.size() && (delim.find(str[end]) == string::npos)) {
      end++; // skip to end of word
    }
    if (end-start != 0) {  // just ignore zero-length strings.
      parts.push_back(string(str, start, end-start));
    }
  }
}

bool getMultiFASTA(vmat_string& fasta, const std::string& filename, bool toUpperCase, bool omitLowerCase){
    char buffer[80];
    ifstream ifs(filename.c_str());
    assert(ifs.fail() == false && "cannot open the file.");

    int index=0;
    fasta.clear();
    fasta.push_back( vec_string(2,"") );

    while( !ifs.eof() )
    {
        // --- name ---
        char c;
        ifs.get(c);
        assert(c == '>');
        std::getline(ifs, fasta[index][0]);
        //fasta[index][0].append((char *)(buffer + 1));
        // --- sequence ---
        fasta[index][1] = "";
        while( !ifs.eof() ){
            char c;
            if(!ifs.get(c))goto MF_OUT;
            if(c == '>'){
                //次の配列データだったとき
                ifs.putback(c);
                index++;
                fasta.push_back( vec_string(2,"") );
                break;
            }
            ifs.putback(c);
            if(!ifs.get(buffer, 80, '>'))goto MF_OUT;
            fasta[index][1].append(buffer);
        }
    }
MF_OUT:
    ifs.close();

    // omit spaces
    for(int i=0; i<(int)fasta.size(); i++)
    {
        string ts="";
        string::iterator s=fasta[i][1].begin();
        while(s!=fasta[i][1].end())
        {
            if(*s != '\r' && *s != '\t' && *s != '\n' && *s != '\f')ts += *s;
            s++;
        }
        fasta[i][1] = ts;
    }

    // to uppercase
    if(toUpperCase)
    {
        for(int j=0; j<(int)fasta.size(); j++)
        {
            for(int i=0; i<(int)fasta.at(j).at(1).size(); i++)
            {
                fasta.at(j).at(1).at(i) = toupper(fasta.at(j).at(1).at(i));
            }
        }
    }

    // omit lowercase region
    if(omitLowerCase)
    {
        // omit spaces
        for(int i=0; i<(int)fasta.size(); i++)
        {
            string ts="";
            string::iterator s=fasta[i][1].begin();
            while(s!=fasta[i][1].end())
            {
                if(! islower(*s))ts += *s;
                s++;
            }
            fasta[i][1] = ts;
        }
    }

    // check for each symbols
    for(int i=0; i<(int)fasta.size(); i++)
    {
        string::iterator s=fasta[i][1].begin();
        while(s!=fasta[i][1].end())
        {
            if(ScoreMatrix::RESIDUE_INDEX[toupper(*s)] == 127)
            {
                string errStr = string(1, (char)toupper(*s)) +
                        " in the sequence " + IS(i+1) +
                        " may be incorrect.";
                putError(errStr);
            }
            s++;
        }
    }


    return true;
}

bool getMultiFASTA(vec_string& names, vec_string& seqs, const std::string& filename, bool toUpperCase, bool omitLowerCase){
    vmat_string fasta;
    getMultiFASTA(fasta, filename, toUpperCase, omitLowerCase);
    names.clear();
    seqs.clear();
    for(vmat_string::iterator i=fasta.begin(); i!=fasta.end(); ++i)
    {
        names.push_back(i->at(0));
        seqs.push_back(i->at(1));
    }

    return true;
}

// flag[i] == true な場所を削除した文字列を返す
std::string* genOmittedString(const std::string& s, const std::vector<bool>& flag) {
    string* newStr = new string();

    //TODO: 一文字ずつ追加で遅いコードをなんとかする
    for(size_t b=0; b<flag.size(); b++)
    {
        if(flag[b] == false) newStr->append(1, s[b]);
    }

    return newStr;
}
