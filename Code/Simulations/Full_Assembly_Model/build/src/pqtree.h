// $Id: pqtree.h 1013 2007-10-08 17:56:50Z cvsrep $
/***************************************************************************
 *   Copyright (C) 2004 by Virginia Tech                                   *
 *   ggrothau@vt.edu                                                       *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/**
PQ-Tree class based on the paper entitled "Testing for the Consecutive Onces Property, Interval Graphs, and Graph Planarity Using PQ-Tree Algorithms" by Kellog S. Booth and George S. Lueker in the Journal of Computer and System Sciences 13, 335-379 (1976)

@author Virginia Tech
*/

#include<list>
#include<set>
#include<vector>
#include<queue>
#include<map>
#include "pqnode.h"
#include "setTemplates.h"
using namespace std;

#ifndef PQTREE_H
#define PQTREE_H

//variable to initiate debugging output
//simply leave this false
bool debugOut=false;

template <class T> class pqtree
{
	private:
	
	
	
	//the root node of the pqtree
	pqnode<T> *root;
	
	//The number of blocks of blocked nodes during the 1st pass
	int block_count;

	//The number of blocked nodes during the 1st pass
	int blocked_nodes;

	//A variable (0 or 1) which is a count of the number of virtual nodes which are
	//imagined to be in the queue during the bubbling up
	int off_the_top;

	//keeps track of all reductions performed on this tree in order
	list<set<T> > reductions;

	//keeps a pointer to the leaf containing T
	//this map actually increases the time complexity of the algorithm
	//to fix, you can create an array of items so that each item hashes
	//to its leaf address in constant time, but this is a tradeoff to
	//conserve space
	map<T,pqnode<T>*> leafAddress;

	//a reference to a pseudonode that cannot be reached through the root
	//of the tree.  The pseudonode is a temporary node designed to handle
	//a special case in the first bubbling up pass
	//it only exists during the scope of the reduce operation
	pqnode<T>* pseudonode;

	//true if a non-safe reduce has failed, tree is useless
	bool invalid;

	//loops through the consecutive blocked siblings of an unblocked 
	//node recursively unblocking the siblings
	int unblockSiblings(pqnode<T>* X,pqnode<T>* parent,pqnode<T>* last);

	/* 
	 * all the different templates for matching a reduce
	 * are below.  The template has a letter describing
	 * which type of node it refers to and a number
	 * /indicating the index of the template for that letter
         * These are the same indeces in the Booth & Lueker paper
	 * The return value indicates whether or not the pattern
	 * accurately matches the template
	 */

	bool template_L1(pqnode<T>*);
	bool template_Q1(pqnode<T>*,bool);
	bool template_Q2(pqnode<T>*,bool);
	bool template_Q3(pqnode<T>*);
	bool template_P1(pqnode<T>*,bool);
	bool template_P2(pqnode<T>*);
	bool template_P3(pqnode<T>*);
	bool template_P4(pqnode<T>*);
	bool template_P5(pqnode<T>*);
	bool template_P6(pqnode<T>*);
	
	//This procedure is the first pass of the Booth&Leuker PQTree algorithm
	//It processes the pertinent subtree of the PQ-Tree to determine the mark
	//of every node in that subtree
	//the pseudonode, if created, is returned so that it can be deleted at
	//the end of the reduce step
	bool bubble(set<T> S);

	bool reduceStep(set<T> S);
				
	public:	
	
	//copy constructor
	pqtree<T>(const pqtree<T>& toCopy);
	//default constructor - constructs a tree using a set
	//only reductions using elements of that set will succeed
	pqtree(set<T> S);
	
	//default destructor
	~pqtree();
 
	//assignment operator
	const pqtree<T>& operator=(const pqtree& toCopy);	
	
	//mostly for debugging purposes, prints the tree to standard out
	void print();
	
	//cleans up pointer mess caused by having a pseudonode
	void cleanPseudo();

	//reduces the tree but protects if from becoming invalid
	//if the reduction fails, takes more time
	bool safeReduce(set<T>);
	bool safeReduceAll(list<set<T> >);

	//reduces the tree - tree can become invalid, making all further
	//reductions fail
	bool reduce(set<T> S);
	bool reduceAll(list<set<T> > L);
	
	//returns 1 possible frontier, or ordering preserving the reductions
	list<T> frontier();
	
	//returns a frontier not including leaves that were not a part of any
	//reduction
	list<T> reducedFrontier();
	
	//returns the reductions that have been performed on this tree
	list<set<T> > getReductions();
	
	//returns the set of all elements on which a reduction was performed with
	set<T> getContained();
};

#endif

//copy constructor
template <class T>
pqtree<T>::pqtree(const pqtree<T>& toCopy)
{
	root=new pqnode<T>(*(toCopy.root));
	block_count=toCopy.block_count;
	blocked_nodes=toCopy.blocked_nodes;
	off_the_top=toCopy.off_the_top;
	invalid=toCopy.invalid;
	reductions=toCopy.reductions;
	pseudonode=NULL;
	leafAddress.clear();
	root->findLeaves(leafAddress);	
}
	
//assignment operator
template <class T>
const pqtree<T>& pqtree<T>::operator=(const pqtree& toCopy)
{
	//make sure we aren't copying ourself
	if(*toCopy==this)
		return *this;
	
	root=new pqnode<T>(*(toCopy.root));
	block_count=toCopy.block_count;
	blocked_nodes=toCopy.blocked_nodes;
	off_the_top=toCopy.off_the_top;
	reductions=toCopy.reductions;
	leafAddress.clear();
	root->findLeaves(leafAddress);	
	
	return this;
}
	
//loops through the consecutive blocked siblings of an unblocked 
//node recursively unblocking the siblings
template <class T>
int pqtree<T>::unblockSiblings(pqnode<T>* X,pqnode<T>* parent,pqnode<T>* last)
{
	int unblockedcnt=0;
	if(X->mark==blocked)
	{
		//unblock the current sibling
		unblockedcnt++;
		X->mark=unblocked;

		X->parent=parent;
		X->check_parent();

		//unblock adjacent siblings recursively
		for(typename set<pqnode<T> *>::iterator i=X->immediate_siblings.begin();i!=X->immediate_siblings.end();i++)
			if(*i!=last)
				unblockedcnt+=unblockSiblings(*i,parent,X);
	}
	return unblockedcnt;
}

//all the different templates for matching a reduce
//are below.  The template has a letter describing
//which type of node it refers to and a number
//indicating the index of the template for that letter
//These are the same indeces in the Booth & Lueker paper

template <class T>
bool pqtree<T>::template_L1(pqnode<T>* X)
{
	//check against the pattern
	if(X->type!=leaf)
		return false;
	
	if(debugOut)
		cerr<<"L1 is running"<<endl;	
	
	X->label=full;
	if(X->parent!=NULL)
		X->parent->full_children.insert(X);
	
	return true;
}
template <class T>
bool pqtree<T>::template_Q1(pqnode<T>* X,bool isRoot)
{	
	//check against pattern
	if(X->type!=qnode)
		return false;

	//find an endmost child
	pqnode<T>* cur=X->endmost_children[0];
	pqnode<T>* last=NULL;	
	//walk the pattern to find any empty/partial children
	while(cur!=NULL)
	{
		if(cur->label!=full)
			return false;
		pqnode<T>* next=cur->qNextChild(last);
		last=cur;
		cur=next;
	}
	
	if(debugOut)
		cerr<<"Q1 is running"<<endl;
	
	X->label=full;		
	X->check_parent();
	if(X->parent!=NULL && !isRoot)
		X->parent->full_children.insert(X);
	return true;
}
	
template <class T>
bool pqtree<T>::template_Q2(pqnode<T>* X,bool isRoot)
{
	//check against the pattern
	if(X->type!=qnode)
		return false;
	if(X->pseudonode)
		return false;
	if(X->partial_children.size()>1)
		return false;
		
	if(X->full_children.size()>0)
	{
		int numfullside=0;
		pqnode<T>* Y;
		pqnode<T>* last=NULL;
		//let Y be the unique element in X's endmost_children that is full
		for(int i=0;i<2;i++)
			if(X->endmost_children[i]->label==full)
			{
				Y=X->endmost_children[i];
				numfullside++;
			}
		//if there are two full endmost_children, our pattern does not match
		if(numfullside!=1)
			return false;
	
		//check that all the full children are consecutive and next to
		//the partial child
		for(int i=0;i<X->full_children.size();i++)
		{
			if(Y->label!=full)
				return false;
			pqnode<T>* next=Y->qNextChild(last);
			last=Y;
			Y=next;
		}
		if(Y->label!=partial && X->partial_children.size()==1)
			return false;
	}
	else
	{
		//or that the partial child is on the side
		int numpartialside=0;
		for(int i=0;i<2;i++)
			if(X->endmost_children[i]->label==partial)
					numpartialside++;
		if(numpartialside==0)
			return false;
	}
	
	if(debugOut)
		cerr<<"Q2 is running"<<endl;

	//if there are no partial children, we are already done
	//otherwise, we need to merge the partial child into X
	if(X->partial_children.size()>0)
	{
		//let Y be the unique partial child of X
		pqnode<T>* Y = *(X->partial_children.begin());
		//let FC be the unique endmost full child of Y
		pqnode<T>* FC;
		//let EC be the unique endmost empty child of Y
		pqnode<T>* EC;
		//find EC and FC
		for(int i=0;i<2;i++)
		{
			if(Y->endmost_children[i]->label==full)	
				FC=Y->endmost_children[i];
			if(Y->endmost_children[i]->label==empty)	
				EC=Y->endmost_children[i];
		}
		
		//let FS be Y's full immediate sibling if Y has one
		pqnode<T>* FS = NULL;
		//let ES be Y's empty immediate sibling if Y has one
		pqnode<T>* ES = NULL;
		//find ES and FS
		for(typename set<pqnode<T> *>::iterator i=Y->immediate_siblings.begin();i!=Y->immediate_siblings.end();i++)
		{
			if((*i)->label==full)
				FS=*i;
			if((*i)->label==empty)
				ES=*i;
		}
		
		//if Y has a full immediate sibling
		if(FS!=NULL)
		{
			FS->immediate_siblings.erase(Y);
			FS->immediate_siblings.insert(FC);
			FC->immediate_siblings.insert(FS);
		}
		else
		{
			if(X->endmost_children[0]==Y)
				X->endmost_children[0]=FC;
			// Since Y has at least two children, X must
			// have at least two children, too, so we
			// should not attach FC on both ends.  Hence
			// the "else" here:
			else if(X->endmost_children[1]==Y)
				X->endmost_children[1]=FC;
			FC->parent=X;
			FC->check_parent();
		}
		// now the FC end of Y's children is spliced into X
		// child list.  The other end is hanging in the air,
		// so is the point in X's child list where Y was taken out.

		//if Y has an empty immediate sibling
		if(ES!=NULL)
		{
			ES->immediate_siblings.erase(Y);
			ES->immediate_siblings.insert(EC);
			EC->immediate_siblings.insert(ES);
		}
		else
		{
			if(X->endmost_children[0]==Y)
				X->endmost_children[0]=EC;
			// Since Y has at least two children, X must
			// have at least two children, too, so we
			// should not attach FS on both ends.  Hence
			// the "else" here:
			else if(X->endmost_children[1]==Y)
				X->endmost_children[1]=EC;
			EC->parent=X;
			EC->check_parent();
		}
	
		//we dont need Y anymore, but we dont want the recursive
		//destructor to kill the children either	
		Y->endmost_children[0]=NULL;	
		Y->endmost_children[0]=NULL;	
		delete Y;
	}
	X->label=partial;	
	if(X->parent!=NULL && !isRoot)
		X->parent->partial_children.insert(X);
	return true;
}
template <class T>
bool pqtree<T>::template_Q3(pqnode<T>* X)
{	
	//check against the pattern
	if(X->type!=qnode)
		return false;
	if(X->partial_children.size()>2)
		return false;

	//handle a special case
	if(X->pseudonode && X->endmost_children[1]==NULL)
		return true;

	bool start=false;
	bool end=false;
	pqnode<T>* cur=X->endmost_children[0];
	pqnode<T>* last=NULL;
	//after this consecutive walk, we know that we have
	//a valid layout

	bool lastIter=false;
	while(!lastIter)
	{
		if(debugOut)
			cerr<<cur<<endl;
		if(cur==X->endmost_children[1])
			lastIter=true;
		if(cur->label==full)
		{
			if(end)
				return false;
			if(!start)
				start=true;
		}
		else if(cur->label==empty)
		{
			if(start)
				end=true;
		}
		else if(cur->label==partial)
		{
			if(!start)
				start=true;
			else if(!end)
				end=true;
			else if(end)
				return false;
		}
		if(!lastIter)
		{
			pqnode<T>* next=cur->qNextChild(last);
			last=cur;
			cur=next;
		}
	}
	
	if(debugOut)
		cerr<<"Q3 is running"<<endl;
	
	//start with the first partial child
	for(typename set<pqnode<T> *>::iterator j=X->partial_children.begin();j!=X->partial_children.end();j++)
	{
		pqnode<T>* PC = *(j);
		pqnode<T>* EC,*FC;	//empty and full child of this partial child
		//find those children
		if(PC->endmost_children[0]->label==full)
		{
			FC=PC->endmost_children[0];
			EC=PC->endmost_children[1];
		}
		else
		{
			EC=PC->endmost_children[0];
			FC=PC->endmost_children[1];
		}
		
		pqnode<T>* CS;
		//handle the sides where the child has immediate siblings
		for(typename set<pqnode<T> *>::iterator i=PC->immediate_siblings.begin();i!=PC->immediate_siblings.end();i++)
		{
			if((*i)->label==empty)
			{
				(*i)->immediate_siblings.erase(PC);
				(*i)->immediate_siblings.insert(EC);
				EC->immediate_siblings.insert(*i);
				CS=FC;
			}
			else	//either full or partial, we dont care which
			{
				(*i)->immediate_siblings.erase(PC);
				(*i)->immediate_siblings.insert(FC);
				FC->immediate_siblings.insert(*i);
				CS=EC;
			}
		}
		if(X->pseudonode)
			CS=EC;
		//handle the case where the child has only one immediate sibling
		if(PC->immediate_siblings.size()==1 || X->pseudonode)
		{
			CS->parent=X;
			CS->check_parent();
			if(X->endmost_children[0]==PC)
				X->endmost_children[0]=CS;
			if(X->endmost_children[1]==PC)
				X->endmost_children[1]=CS;
		}
		
		//we want to delete PC, but not it's children
		PC->endmost_children[0]=NULL;
		PC->endmost_children[1]=NULL;
		delete PC;
		PC=NULL;

	}
	return true;
}
/*A note here.  An error in the Booth and Leuker
Algorithm fails to consider the case where a P-node
is full, is the pertinent root, and is not an endmost
child of a q-node.  In this case, we need to know that
the P-node is a pertinent root and not try to update its
parent whose pointer is possibly invalid*/
template <class T>
bool pqtree<T>::template_P1(pqnode<T>* X,bool isRoot)
{
		//check against the pattern
		if(X->type!=pnode)
			return false;
		if(X->full_children.size()!=X->child_count())
			return false;
		if(debugOut)
			cerr<<"P1 is running"<<endl;

		X->label=full;
		if(!isRoot)	//make sure we aren't root first
			X->parent->full_children.insert(X);
		
		return true;
	}
template <class T>
bool pqtree<T>::template_P2(pqnode<T>* X)
{
	//check against the pattern
	if(X->type!=pnode)
		return false;
	if(X->partial_children.size()!=0)
		return false;
		
	if(debugOut)
			cerr<<"P2 is running"<<endl;

	//move X's full children into their own P-node
	if(X->full_children.size()>=2)
	{
		pqnode<T>* Y=new pqnode<T>;
		Y->parent=X;
		Y->check_parent();
		Y->type=pnode;
		for(typename set<pqnode<T> *>::iterator i=X->full_children.begin();i!=X->full_children.end();i++)
		{
			X->circular_link.remove(*i);
			Y->circular_link.push_back(*i);
			(*i)->parent=Y;
			(*i)->check_parent();
		}
		X->circular_link.push_back(Y);
	}
	//mark the root partial
	X->label=partial;
		
	return true;
}
template <class T>
bool pqtree<T>::template_P3(pqnode<T>* X)
{
	//check against the pattern
	if(X->type!=pnode)
		return false;
	if(X->partial_children.size()!=0)
		return false;
		
	if(debugOut)
		cerr<<"P3 is running"<<endl;

	pqnode<T>* theParent=X->parent;

	//create new Q-node Y
	pqnode<T>* Y = new pqnode<T>;
	Y->parent=theParent;
	Y->check_parent();
	Y->type=qnode;

	//switch Y for X as parent's children
	theParent->partial_children.insert(Y);
	theParent->partial_children.erase(X);
	if(theParent->type==pnode)
	{
		theParent->circular_link.remove(X);
		theParent->circular_link.push_back(Y);
	}
	else
		X->swapQ(Y);

	//set up the p_node on the full child side
	pqnode<T>* FC,*EC;
	if(X->full_children.size()==1)
	{
		FC=*(X->full_children.begin());
		X->circular_link.remove(FC);
	}
	else
	{
		FC = new pqnode<T>;	//FC = full child
		FC->label=full;
		FC->type=pnode;
		for(typename set<pqnode<T> *>::iterator i=X->full_children.begin();i!=X->full_children.end();i++)
		{
			X->circular_link.remove(*i);
			FC->circular_link.push_back(*i);
			(*i)->parent=FC;
			(*i)->check_parent();
		}
	}
	FC->parent=Y;
	FC->check_parent();
	Y->endmost_children[0]=FC;
	Y->full_children.insert(FC);

	//now set up the p-node on the empty child side
	if(X->circular_link.size()==1)
	{
		EC=*(X->circular_link.begin());

		//we want to delete X, but not it's children
		X->circular_link.clear();
		delete X;
	}
	else
		EC=X;
	EC->parent=Y;
	EC->check_parent();
	EC->label=empty;
	Y->endmost_children[1]=EC;

	//update the immediate siblings links
	EC->immediate_siblings.clear();
	EC->immediate_siblings.insert(FC);
	FC->immediate_siblings.clear();
	FC->immediate_siblings.insert(EC);
	
	Y->label=partial;
	
	return true;
}
template <class T>
bool pqtree<T>::template_P4(pqnode<T>* X)
{
	//check against the pattern
	if(X->type!=pnode)
		return false;
	if(X->partial_children.size()!=1)
		return false;


	
	//Y is the partial Q-node
	pqnode<T>* Y=*X->partial_children.begin();
	pqnode<T>* EC;
	pqnode<T>* FC;
	pqnode<T>* ES;  //empty/full child/sibling of Y

	//find the empty and full endmost child of Y
	if(Y->endmost_children[0]->label==empty)
	{
		EC=Y->endmost_children[0];
		FC=Y->endmost_children[1];
	}
	else
	{
		EC=Y->endmost_children[1];
		FC=Y->endmost_children[0];
	}
	//check that we are indeed matching the pattern
	if(EC->label!=empty || FC->label==empty)
		return false;
	
	if(debugOut)
		cerr<<"P4 is running"<<endl;
	
	//if Y has an empty sibling, set ES to be an empty sibling of Y
	for(typename list<pqnode<T> *>::iterator i=X->circular_link.begin();i!=X->circular_link.end();i++)
		if((*i)->label==empty)
		{
			if(debugOut)
				cerr<<"Y has an empty sibling"<<endl;
			ES=*i;
			break;
		}
	
	//move the full children of X to be children of Y
	if(X->full_children.size()>0)
	{
		pqnode<T> *ZF;
		//only 1 full child
		if(X->full_children.size()==1)
		{
			if(debugOut)
				cerr<<"Only 1 full child"<<endl;
			
			ZF=*(X->full_children.begin());
			X->circular_link.remove(ZF);
		}
		//multiple full children - must be placed in a P-node
		else
		{
			//create ZF to be a new p-node
			ZF=new pqnode<T>;
			ZF->label=full;
			ZF->type=pnode;

			//for all the full nodes in X, set their parent to be ZF
			for(typename set<pqnode<T> *>::iterator W=X->full_children.begin();W!=X->full_children.end();W++)
			{
				(*W)->parent=ZF;
				(*W)->check_parent();
				X->circular_link.remove(*W);
				ZF->circular_link.push_back(*W);
			}
		}
		//more updates
		ZF->parent=Y;
		ZF->check_parent();
		FC->immediate_siblings.insert(ZF);			
		ZF->immediate_siblings.insert(FC);
		if(Y->endmost_children[0]==FC)
		  	Y->endmost_children[0]=ZF;
                else if(Y->endmost_children[1]==FC)
			Y->endmost_children[1]=ZF;
		Y->full_children.insert(ZF);
	}

	//if X now only has one child, get rid of X
	if(X->circular_link.size()==1)
	{
		if(debugOut)
			cerr<<"X only has one child"<<endl;
		pqnode<T>* theParent = X->parent;
		Y->parent = X->parent;
		
		if(theParent!=NULL)	//parent is root of tree
		{

			//update parent to handle the switch
			if(X->immediate_siblings.empty())	//parent is a p-node
			{
				theParent->circular_link.remove(X);	
				theParent->circular_link.push_back(Y);
			}
			else	//parent is a Q-node
			{
				//update the immediate siblings list by removing X and adding Y
				for(typename set<pqnode<T> *>::iterator i=X->immediate_siblings.begin();i!=X->immediate_siblings.end();i++)
				{
					(*i)->immediate_siblings.erase(X);
					(*i)->immediate_siblings.insert(Y);
					Y->immediate_siblings.insert(*i);
				}
				if(X->immediate_siblings.size()==1)
				{
					if(theParent->endmost_children[0]==X)
						theParent->endmost_children[0]=Y;
					if(theParent->endmost_children[1]==X)
						theParent->endmost_children[1]=Y;
				}
			}
		}
		else
			root=Y;
	}
	
	return true;
}
int runs=0;
template <class T>
bool pqtree<T>::template_P5(pqnode<T>* X)
{	
	//check against the pattern
	if(X->type!=pnode)
		return false;
	if(X->partial_children.size()!=1)
		return false;
		
	
	//Y is the partial Q-node
	pqnode<T>* Y=*X->partial_children.begin();
	pqnode<T>* EC=NULL;
	pqnode<T>* FC=NULL;
	pqnode<T>* ES=NULL;  //empty/full child/sibling of Y

	//find the empty and full endmost child of Y
	if(Y->endmost_children[0]->label==empty)
	{
		EC=Y->endmost_children[0];
		FC=Y->endmost_children[1];
	}
	else
	{
		EC=Y->endmost_children[1];
		FC=Y->endmost_children[0];
	}
	//check that we are indeed matching the pattern
	if(EC->label!=empty || FC->label==empty)
		return false;
		
	if(debugOut)
		cerr<<"P5 is running"<<endl;
		
	//if Y has an empty sibling, set ES to be an empty sibling of Y
	for(typename list<pqnode<T> *>::iterator i=X->circular_link.begin();i!=X->circular_link.end();i++)
		if((*i)->label==empty)
		{
			ES=*i;
			break;
		}

	//Y will be the root of the pertinent subtree after the replacement
	pqnode<T>* theParent = X->parent;
	Y->parent = X->parent;
	Y->pertinent_leaf_count = X->pertinent_leaf_count;
	Y->label=partial;
	//add Y to it's parent's list of partial children
	theParent->partial_children.insert(Y);
	//remove Y from X's list of children
	X->circular_link.remove(Y);
	X->partial_children.erase(Y);
	
	//update parent to handle the switch
	if(X->immediate_siblings.size()==0)	//parent is a P-node
	{
		theParent->circular_link.remove(X);
		theParent->circular_link.push_back(Y);
	}
	else	//parent is a Q-node
	{
		//update the immediate siblings list by removing X and adding Y
		for(typename set<pqnode<T> *>::iterator i=X->immediate_siblings.begin(); i!=X->immediate_siblings.end();i++)
		{
			(*i)->immediate_siblings.erase(X);
			(*i)->immediate_siblings.insert(Y);
			Y->immediate_siblings.insert(*i);
		}
		if(theParent->endmost_children[0]==X)
			theParent->endmost_children[0]=Y;
		if(theParent->endmost_children[1]==X)
			theParent->endmost_children[1]=Y;
	}
	
	//move the full children of X to be children of Y
	if(X->full_children.size()>0)
	{
		pqnode<T> *ZF=NULL;
		//only 1 full child
		if(X->full_children.size()==1)
		{
			ZF=*(X->full_children.begin());
			X->circular_link.remove(ZF);
		}
		//multiple full children - must be placed in a P-node
		else
		{
			//create ZF to be a new p-node
			ZF=new pqnode<T>;
			ZF->label=full;
			ZF->type=pnode;

			//for all the full nodes in X, set their parent to be ZF
			for(typename set<pqnode<T> *>::iterator W=X->full_children.begin();W!=X->full_children.end();W++)
			{
				(*W)->parent=ZF;
				(*W)->check_parent();
				X->circular_link.remove(*W);
				ZF->circular_link.push_back(*W);
			}
		}
		X->full_children.clear();
		
		//more updates
		ZF->parent=Y;
		ZF->check_parent();
		FC->immediate_siblings.insert(ZF);			
		ZF->immediate_siblings.insert(FC);
		if(Y->endmost_children[0]==FC)
                	Y->endmost_children[0]=ZF;
                else if(Y->endmost_children[1]==FC)
                        Y->endmost_children[1]=ZF;			
	}

	//if X still has some empty children, insert them	
	if(X->child_count()>0)
	{
		pqnode<T> *ZE=NULL;
		if(X->child_count()==1)
			ZE=ES;
		else
		{
			ZE=X;
			ZE->label=empty;
			ZE->immediate_siblings.clear();
		}
		ZE->parent=Y;
		ZE->check_parent();
		EC->immediate_siblings.insert(ZE);
		ZE->immediate_siblings.insert(EC);
		if(Y->endmost_children[0]==EC)
			Y->endmost_children[0]=ZE;
		if(Y->endmost_children[1]==EC)
			Y->endmost_children[1]=ZE;
	}
	if(X->child_count()<2)
	{
		//we want to delete X, but not it's children
		X->circular_link.clear();
		delete X;		
	}
	

	return true;
}

template <class T>
bool pqtree<T>::template_P6(pqnode<T>* X)
{
	//check against the pattern
	if(X->type!=pnode)
		return false;
	if(X->partial_children.size()!=2)
		return false;

	
	//Y is the first partial Q-node from which we shall build
	pqnode<T>* Y=*X->partial_children.begin();
	pqnode<T>* Z=*(++(X->partial_children.begin()));
	pqnode<T>* EC=NULL;
	pqnode<T>* FC=NULL;  //empty/full child/sibling of Y

	pqnode<T> *ZF=NULL;	//the child of Y created to hold full X's children

	//find the empty and full endmost child of Y
	if(Y->endmost_children[0]->label==empty)
	{
		EC=Y->endmost_children[0];
		FC=Y->endmost_children[1];
	}
	else
	{
		EC=Y->endmost_children[1];
		FC=Y->endmost_children[0];
	}
	//check that we are indeed matching the pattern
	if(EC->label!=empty || FC->label==empty)
		return false;
		
	if(debugOut)
		cerr<<"P6 is running"<<endl;

	//move the full children of X to be children of Y
	if(X->full_children.size()>0)
	{
		//only 1 full child
		if(X->full_children.size()==1)
		{
			ZF=*(X->full_children.begin());
			X->circular_link.remove(ZF);
		}
		//multiple full children - must be placed in a P-node
		else
		{
			//create ZF to be a new p-node
			ZF=new pqnode<T>;
			ZF->label=full;
			ZF->type=pnode;

			//for all the full nodes in X, set their parent to be ZF
			for(typename set<pqnode<T> *>::iterator W=X->full_children.begin();W!=X->full_children.end();W++)
			{
				(*W)->parent=ZF;
				(*W)->check_parent();
				X->circular_link.remove(*W);
				ZF->circular_link.push_back(*W);
			}
		}
		//more updates
		ZF->parent=Y;
		ZF->check_parent();
		FC->immediate_siblings.insert(ZF);			
		ZF->immediate_siblings.insert(FC);
		if(Y->endmost_children[0]==FC)
	       		Y->endmost_children[0]=ZF;
		if(Y->endmost_children[1]==FC)
			Y->endmost_children[1]=ZF;			
		
		//now, incorporate the other partial child

		//find the empty and full endmost child of Z
		if(Z->endmost_children[0]->label==empty)
		{
			EC=Z->endmost_children[0];
			FC=Z->endmost_children[1];
		}
		else
		{
			EC=Z->endmost_children[1];
			FC=Z->endmost_children[0];
		}
		//check that we are indeed matching the pattern
		if(EC->label!=empty || FC->label==empty)
			return false;

		//connect the children of the two partials together
		ZF->immediate_siblings.insert(FC);
		FC->immediate_siblings.insert(ZF);

		//adjust the parent Y
		if(Y->endmost_children[0]==ZF)
			Y->endmost_children[0]=EC;
		if(Y->endmost_children[1]==ZF)
			Y->endmost_children[1]=EC;
		FC->parent=Y;
		FC->check_parent();
		EC->parent=Y;
		EC->check_parent();
	}
	else	//no full children
	{
		pqnode<T> *ZFC=NULL;	//Z's full child

		//figure out which sides of Y and Z to connect together
		if(Z->endmost_children[0]->label==full)
		{
			ZFC=Z->endmost_children[0];
			if(Y->endmost_children[0]==FC)
				Y->endmost_children[0]=Z->endmost_children[1];
			if(Y->endmost_children[1]==FC)
				Y->endmost_children[1]=Z->endmost_children[1];
		}
		else
		{
			ZFC=Z->endmost_children[1];
			if(Y->endmost_children[0]==FC)
				Y->endmost_children[0]=Z->endmost_children[0];
			if(Y->endmost_children[1]==FC)
				Y->endmost_children[1]=Z->endmost_children[0];
		}
		Y->endmost_children[0]->parent=Y;
		Y->endmost_children[0]->check_parent();
		Y->endmost_children[1]->parent=Y;
		Y->endmost_children[1]->check_parent();

		FC->immediate_siblings.insert(ZFC);
		ZFC->immediate_siblings.insert(FC);
	}
		
	
	//adjust the root X
	//we dont need Z any more
	X->circular_link.remove(Z);
	Z->endmost_children[0]=NULL;
	Z->endmost_children[1]=NULL;
	delete Z;

	//if X now only has one child, get rid of X
	if(X->circular_link.size()==1)
	{
	
		pqnode<T>* theParent = X->parent;
		Y->parent = X->parent;
		Y->pertinent_leaf_count = X->pertinent_leaf_count;
		Y->label=partial;
	
		if(theParent!=NULL)	//parent is not root of tree
		{	
			//add Y to it's parent's list of partial children
			theParent->partial_children.insert(Y);
	
			//update parent to handle the switch
			if(theParent->type==pnode)
			{
				theParent->circular_link.remove(X);	
				theParent->circular_link.push_back(Y);
			}
			else	//parent is a Q-node
			{
				//update the immediate siblings list by removing X and adding Y
				for(typename set<pqnode<T> *>::iterator i=X->immediate_siblings.begin();i!=X->immediate_siblings.end();i++)
				{
					(*i)->immediate_siblings.erase(X);
					(*i)->immediate_siblings.insert(Y);
				}
				if(theParent->endmost_children[0]==X)
					theParent->endmost_children[0]=Y;
				if(theParent->endmost_children[1]==X)
					theParent->endmost_children[1]=Y;
			}
		}
		else	//this is the root of our tree, act accordingly
		{
			root=Y;
			Y->parent=NULL;
			
			//delete X, but not it's children
			X->circular_link.clear();
			delete X;
		}	
	}
	return true;
}

//This procedure is the first pass of the Booth&Leuker PQTree algorithm
//It processes the pertinent subtree of the PQ-Tree to determine the mark
//of every node in that subtree
//the pseudonode, if created, is returned so that it can be deleted at
//the end of the reduce step
template <class T>
bool pqtree<T>::bubble(set<T> S)
{
	if(debugOut) {
	cerr<<"Running bubble on set ";
	for(typename set<T>::iterator i=S.begin();i!=S.end();i++)
		cerr<<*i<<",";
	cerr<<endl;}
	//initialize variables
	queue<pqnode<T>*> q;
	block_count=0;
	blocked_nodes=0;
	off_the_top=0;
	
	set<pqnode<T>* > blocked_list;	//stores blocked guys that arent unblocked

	//insert the set's leaves into the queue
	for(typename set<T>::iterator i=S.begin();i!=S.end();i++)
	{
		pqnode<T> *temp=leafAddress[*i];
		if(temp==NULL)	//make sure we have this in our leaves already
			return false;
		q.push(temp);
	}

	while(q.size()+block_count+off_the_top>1)
	{
		//check to see if there are still guys in the queue
		if(q.empty())
			return false;
		
		//remove X from the front of the queue
		pqnode<T>* X = q.front();
		q.pop();
		

		//mark X as blocked
		X->mark=blocked;
		
		//get the set of blocked and unblocked siblings
		set<pqnode<T> *> US, BS;
		for(typename set<pqnode<T> *>::iterator i=X->immediate_siblings.begin();i!=X->immediate_siblings.end();i++)
		{
			if((*i)->mark==blocked)
				BS.insert(*i);
			if((*i)->mark==unblocked)
				US.insert(*i);
		}
		//we can assign a parent to X if one of its immediate siblings is unblocked
		//or it has 0/1 immediate siblings meaning it is a corner child of a q node
		//or a child of a p node
		if(!US.empty())
		{
			X->parent=(*US.begin())->parent;
			X->check_parent();
			X->mark=unblocked;
		}
		else if(X->immediate_siblings.size()<2)
				X->mark=unblocked;
		
		//if it is unblocked, we can process it
		if(X->mark==unblocked)
		{
			int listSize=0;
			pqnode<T>* Y=X->parent;
			if(!BS.empty())
			{
				X->mark=blocked;	//will be unblocked by recursive unblocksiblings
				listSize=unblockSiblings(X,Y,NULL);
				Y->pertinent_child_count+=listSize-1;
			}
			
			if(Y==NULL)	//currently at root node
				off_the_top=1;
			else
			{
				Y->pertinent_child_count++;
				if(Y->mark==unmarked)
				{
					q.push(Y);
					Y->mark=queued;
				}
			}
			block_count-=BS.size();
			blocked_nodes-=listSize;					
		}
		else
		{
			block_count+=1-BS.size();
			blocked_nodes+=1;
			blocked_list.insert(X);
		}
	}

	if(block_count>1 || (off_the_top==1 && block_count!=0))
		return false;

	
	int correctblockedcount=0;
	for(typename set<pqnode<T>*>::iterator i=blocked_list.begin();i!=blocked_list.end();i++)
			if((*i)->mark==blocked)
				correctblockedcount++;

	//in this case, we have a block that is contained within a Q-node
	//we must assign a psuedonode Z of type Q-node to handle it
	if(block_count==1 && correctblockedcount>1)
	{
		pseudonode=new pqnode<T>;
		if(debugOut)
			cerr<<"pseudonode"<<pseudonode<<endl;
		pseudonode->type=qnode;
		pseudonode->pseudonode=true;
		int side=0;
		pseudonode->pertinent_child_count=0;
		
		//figure out which nodes are still blocked
		//and which of those are the endmost children
		for(typename set<pqnode<T>*>::iterator i=blocked_list.begin();i!=blocked_list.end();i++)
			if((*i)->mark==blocked)
			{
				if(debugOut)
					cerr<<"blocked: "<<*i<<endl;
				pseudonode->pertinent_child_count++;
				pseudonode->pertinent_leaf_count+=(*i)->pertinent_leaf_count;
				int count=0;	//count number of immediate siblings
				int loop=0;
				for(typename set<pqnode<T>*>::iterator j=(*i)->immediate_siblings.begin();j!=(*i)->immediate_siblings.end() && loop<2;j++)
				{
					loop++;

					if((*j)->mark==blocked)
						count++;
					else
					{
						(*i)->immediate_siblings.erase(*j);
						(*j)->immediate_siblings.erase(*i);
						pseudonode->pseudo_neighbors[side]=*j;
					}
				}
				(*i)->parent=pseudonode;
				(*i)->check_parent();
				(*i)->pseudochild=true;
				if(count<2)
					pseudonode->endmost_children[side++]=*i;
			}
	}
	return true;
}//bool bubble(set<T> S)

template <class T>
bool pqtree<T>::reduceStep(set<T> S)
{
	if(debugOut) {
	cerr<<"Running reduceStep on set ";
	for(typename set<T>::iterator i=S.begin();i!=S.end();i++)
		cerr<<*i<<",";
	cerr<<endl;}
	
	//build a queue with all the pertinent leaves in it
	queue<pqnode<T>*> q;
	for(typename set<T>::iterator i=S.begin();i!=S.end();i++)
	{
		pqnode<T>* X = leafAddress[*i];
		if(X==NULL)
			return false;
		X->pertinent_leaf_count=1;
		q.push(X);
	}
	
	while(!q.empty())
	{
		//remove X from the front of the queue
		pqnode<T>* X = q.front();
		q.pop();
		if(debugOut) cerr<<"front of queue "<<X<<endl;

		if(X->pertinent_leaf_count<S.size())	//X is not root
		{
			if(debugOut)	cerr<<"Not Root - "<<X<<endl;
			//update X's parent Y
			pqnode<T>* Y = X->parent;
			Y->pertinent_leaf_count+=X->pertinent_leaf_count;
			Y->pertinent_child_count--;
			//push Y onto the queue if it has no more pertinent children
			if(Y->pertinent_child_count==0)
				q.push(Y);
			
			//testagainst various templates
			if(template_L1(X)) {if(debugOut) cerr<<"Running L1"<<endl;}
			
			else if(template_P1(X,false)) {if(debugOut) cerr<<"Running P1"<<endl;}
			else if(template_P3(X)) {if(debugOut) cerr<<"Running P3"<<endl;}
			else if(template_P5(X)) {if(debugOut) cerr<<"Running P5"<<endl;}
		
			else if(template_Q1(X,false)) {if(debugOut) cerr<<"Running Q1"<<endl;}
			else if(template_Q2(X,false)) {if(debugOut) cerr<<"Running Q2"<<endl;}
			
			else
			{
				cleanPseudo();
				return false;
			}
		}
		else					//X is root
		{	
			if(debugOut)	cerr<<"Root - "<<X<<endl;
			if(template_L1(X)) {if(debugOut) cerr<<"Running L1"<<endl;}
		
			else if(template_P1(X,true)) {if(debugOut) cerr<<"Running P1"<<endl;}
			else if(template_P2(X)) {if(debugOut) cerr<<"Running P2"<<endl;}
			else if(template_P4(X)) {if(debugOut) cerr<<"Running P4"<<endl;}
			else if(template_P6(X)) {if(debugOut) cerr<<"Running P6"<<endl;}
		
			else if(template_Q1(X,true)) {if(debugOut) cerr<<"Running Q1"<<endl;}
			else if(template_Q2(X,true)) {if(debugOut) cerr<<"Running Q2"<<endl;}
			else if(template_Q3(X)) {if(debugOut) cerr<<"Running Q3"<<endl;}
		
			else
			{
				cleanPseudo();
				return false;
			}
		}
	}
	cleanPseudo();
	return true;
}

template <class T>
void pqtree<T>::cleanPseudo()
{
	if(pseudonode!=NULL)
	{
			for(int i=0;i<2;i++)
			{
					pseudonode->endmost_children[i]->immediate_siblings.insert(pseudonode->pseudo_neighbors[i]);
				pseudonode->pseudo_neighbors[i]->immediate_siblings.insert(pseudonode->endmost_children[i]);
			}

			pseudonode->endmost_children[0]=NULL;
			pseudonode->endmost_children[1]=NULL;
			delete pseudonode;
			pseudonode=NULL;
	}
}

//basic constructor from an initial set
template <class T>
pqtree<T>::pqtree(set<T> S)
{
	//set up the root node as a P-Node initially
	root=new pqnode<T>;
	root->type=pnode;	
	invalid=false;
	pseudonode=NULL;
	block_count=0;
	blocked_nodes=0;
	off_the_top=0;
	for(typename set<T>::iterator i=S.begin();i!=S.end();i++)
	{
		pqnode<T> *newNode;
		newNode = new pqnode<T>(*i);
		leafAddress[*i]=newNode;
		newNode->parent = root;
		newNode->type=leaf;
		root->circular_link.push_back(newNode);
	}
}

template <class T>
void pqtree<T>::print()
{
	root->print();
}

//reduces the tree but protects if from becoming invalid
//if the reduction fails, takes more time
template <class T>
bool pqtree<T>::safeReduce(set<T> S)
{
	//using a backup copy to enforce safety
	pqtree<T> toCopy(*this);
	
	if(!reduce(S))
	{
		//reduce failed, so perform a copy
		root=new pqnode<T>(*toCopy.root);
		block_count=toCopy.block_count;
		blocked_nodes=toCopy.blocked_nodes;
		off_the_top=toCopy.off_the_top;
		invalid=toCopy.invalid;
		leafAddress.clear();
		root->findLeaves(leafAddress);	
		return false;
	}
	return true;
}

template <class T>
bool pqtree<T>::safeReduceAll(list<set<T> > L)
{
	//using a backup copy to enforce safety
	pqtree<T> toCopy(*this);
	if(!reduceAll(L))
	{
		//reduce failed, so perform a copy
		root=new pqnode<T>(*toCopy.root);
		block_count=toCopy.block_count;
		blocked_nodes=toCopy.blocked_nodes;		off_the_top=toCopy.off_the_top;
		invalid=toCopy.invalid;
		leafAddress.clear();
		root->findLeaves(leafAddress);	
		return false;
	}
	return true;
}	

	
template <class T>
bool pqtree<T>::reduce(set<T> S)
{
	if(S.size()<2)
	{
		reductions.push_back(S);
		return true;
	}
	if(invalid)
		return false;
	if(!bubble(S))
	{
		invalid=true;
		return false;
	}
	if(!reduceStep(S))
	{
		invalid=true;
		return false;
	}
	
	//we dont need any more pseudonodes
	//at least until the next reduce
	if(pseudonode!=NULL)
	{
		pseudonode->endmost_children[0]=NULL;
		pseudonode->endmost_children[1]=NULL;
		delete pseudonode;
		pseudonode=NULL;
	}

	//reset all the temporary variables for the next round
	root->reset();
	
	//store the reduction set for later lookup
	reductions.push_back(S);

	
	return true;
}

template <class T>
bool pqtree<T>::reduceAll(list<set<T> > L)
{
	if(debugOut)
		cerr<<"Running ReduceAll"<<endl;
	if(invalid)
		return false;

	for(typename list<set<T> >::iterator S=L.begin();S!=L.end();S++)
	{
		if(S->size()<2)
		{
			reductions.push_back(*S);
			continue;
		}
		if(debugOut)
		{
			cerr<<"*"<<endl;
			print();
			cerr<<"%"<<endl;
		}
		if(!bubble(*S))
		{
			invalid=true;
			return false;
		}
		if(!reduceStep(*S))
		{
			invalid=true;
			return false;
		}

		//we dont need any more pseudonodes
		//at least until the next reduce
		if(pseudonode!=NULL)
		{
			pseudonode->endmost_children[0]=NULL;
			pseudonode->endmost_children[1]=NULL;
			delete pseudonode;
			pseudonode=NULL;
		}
		//reset all the temporary variables for the next round
		root->reset();
	
		//store the reduction set for later lookup
		reductions.push_back(*S);
	}	
	return true;

}

template <class T>
list<T> pqtree<T>::frontier()
{
	list<T> out;
	root->findFrontier(out);
	return out;
}

template <class T>
list<T> pqtree<T>::reducedFrontier()
{
	list<T> out,inter;
	root->findFrontier(inter);
	set<T> allContained;
	for(list<set<int> >::iterator j=reductions.begin();j!=reductions.end();j++)
			allContained=setunion(allContained,*j);
	for(list<int>::iterator j=inter.begin();j!=inter.end();j++)
			if(setfind(allContained,*j))
				out.push_back(*j);

	return out;
}

template <class T>
list<set<T> > pqtree<T>::getReductions()
{
	return reductions;
}

template <class T>
set<T> pqtree<T>::getContained()
{
	set<T> out;
	for(typename list<set<T> >::iterator i=reductions.begin();i!=reductions.end();i++)
		out=setunion(out,*i);
	
	return out;
}

//default destructor, just needs to delete the root
//for safety's sake, we delete the pseudonode too
template <class T>
pqtree<T>::~pqtree()
{
	delete root;
	//delete pseudonode;
}
