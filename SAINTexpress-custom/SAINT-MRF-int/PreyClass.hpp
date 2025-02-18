/*
 * baitClass.hpp
 *
 *  Created on: June 18, 2012
 *      Author: hwchoi
 */


#ifndef PREYCLASS_HPP_
#define PREYCLASS_HPP_

#include <iostream>
#include <string>
#include <map>
#include <deque>



using namespace std;



class PreyClass {

private:

	int rowId;
	string preyId;
	double preyLength;
	string preyGeneId;
	double imputedMu;
	


public:
	
	PreyClass();
	void print() const;

	int get_rowId() const;
	string get_preyId() const;
	double get_preyLength() const;
	string get_preyGeneId() const;
	double get_imputedMu() const;

	void set_rowId(int r);
	void set_preyId(string str);
	void set_preyLength(double x);
	void set_preyGeneId(string str);
	void set_imputedMu(double x);

};


#endif /* PREYCLASS_HPP_ */
