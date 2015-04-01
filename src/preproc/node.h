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
 * $Id: node.h 12 2003-11-07 17:58:49Z kdalbey $ 
 */

// node.h for the definition of node class
#ifndef NODE_H
#define NODE_H


class Node {


friend class Element;
//friend class Edge;
friend class Boundary;
friend void get_key(double* , unsigned* , unsigned* );

 public:

  Node(int, double*, int, int); //constructor
  Node();
  void setparameters(int, double*);
  int get_nodeid() {return nodeid;};
  double* get_node_coord(){return node_coord;};
  double* get_coord(){return node_coord;};
  void determine_the_key(unsigned, double*, double*, unsigned*, unsigned*);
  void determine_max_min(double*, double*);
  int get_written_flag(){return written;};
  void set_written_flag(){written=1;};
  unsigned* get_key(){return key;};
  void write_node_data(ofstream*);
  void write_node_data_bin(FILE *);
  void clear_written_flag();
  void put_element_array_loc(int);
  int get_element_array_loc(int);

 private:
  int nodeid;
  double node_coord[3];
  unsigned key[2];
  int written;
  int element_array_loc[2]; // used to find the elements faster

};
#endif
