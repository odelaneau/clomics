/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "mbed_data.h"

void mbed_data::tagModuled(mcluster * curr_node, int curr_mod) {
	int _curr_mod = curr_mod;
	if (curr_node->leaf()) curr_node->idx_module = curr_mod;
	else {
		if (curr_mod < 0 && module_set.count(curr_node->sid) > 0) _curr_mod = curr_node->index;
		tagModuled(curr_node->child1, _curr_mod);
		tagModuled(curr_node->child2, _curr_mod);
	}
}
