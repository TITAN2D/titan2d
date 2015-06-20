/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: element.h 117 2007-06-07 19:14:40Z dkumar $ 
 */

#ifndef ELEMENT_H
#define ELEMENT_H
#include "node.h"
#include "boundary.h"

class ElementPreproc {

	friend void drawit(ElementPreproc);

public:

	ElementPreproc();
	void setparameters(int, NodePreproc*[], int, int*);
	void order_nodes();
	void case1();
	void case2();
	void case3();
	void case4();
	void case5();
	void create_m_node(double*, double*);
	int getid() {
		return elementid;
	}
	;
	//void determine_the_key(unsigned, double*, double*);
	unsigned* pass_key();
	void determine_neighbors(int, ElementPreproc*);
	NodePreproc** get_element_node() {
		return element_nodes;
	}
	;
	Boundary* get_element_bound(int iside, int i) {
		return boundary[iside][i];
	}
	;
	void myproc(int, int, int);
	int get_myprocess() {
		return myprocess;
	}
	;
	int* get_elm_loc() {
		return elm_loc;
	}
	;
	int get_elm_mat() {
		return material;
	}
	;
	void write_node_data(ofstream*);
	void write_element_data(ofstream*);
	void write_element_data_bin(FILE *);
	void reset_written_flag();
	void find_boundary(int, int, Boundary*);
	void set_boundary(Boundary *b);

	//Edge* get_element_edges();
	ElementPreproc* get_neighbors(int);
	void determine_opposite_brother();
	ElementPreproc* get_opposite_brother() {
		return opposite_brother;
	}
	;

private:
	int elementid;
	Boundary *boundary[4][2];
	NodePreproc* element_nodes[9];
	//Edge* element_edges[4];
	ElementPreproc* neighbor[4];
	int material;
	int myprocess;
	int elm_loc[2];
	int which_son;
	ElementPreproc* opposite_brother;

};
#endif
