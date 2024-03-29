// $Id: setTemplates.h 1010 2007-10-08 10:47:19Z cvsrep $
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
 
 /*
  * A group of templates that define commonly used operations on stl sets
  */
 
 
#include<algorithm>
#include<set>
using namespace std;

template <class T> set<T> setintersection(set<T> A, set<T> B)
{
		set<T> out;
		set_intersection(A.begin(),A.end(),B.begin(),B.end(),
			inserter(out,out.begin()) );
		return out;
}
template <class T> set<T> setunion(set<T> A, set<T> B)
{
		set<T> out;
		set_union(A.begin(),A.end(),B.begin(),B.end(),
			inserter(out,out.begin()) );
		return out;
}
template <class T> bool setfind(set<T> A, T B)
{
		set<T> temp;
		temp.insert(B);
		return (setintersection(A,temp).size()>0);
}
template <class T> set<T> setdifference(set<T> A, set<T> B)
{
		set<T> out;
	        set_difference(A.begin(),A.end(),B.begin(),B.end(),
		        inserter(out,out.begin()) );
	        return out;
}

