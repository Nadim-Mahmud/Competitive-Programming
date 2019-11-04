## STL
  * **Multiset**
    * `multiset <int, greater <int> > set1;` empty multiset container greater for descending order
	 * `set1.insert(40);` insert elements in random order
    * `set1.erase(set1.begin(), set1.find(40));` remove all elements up to element with value 40 in set1
	* `int num = set1.erase(50);` remove all elements with value 50 in set2 || num = number of element 50 removed
	* `set1.clear();` clear all elements
    * `set1.find(20);` Returns iterator to 20 if found else returns end()
    * `set1.count(const g) ` Returns the number of matches to element ‘g’ in the multiset.
    * `set1.rbegin()` Returns a reverse iterator pointing to the last element in the multiset container.
    * `multiset <int, greater <int> > :: iterator itr;` iterator is normal like others iterators
	 * `multiset <int> set2(set1.begin(), set1.end());` assigning the elements from set1 to set2	
    * `multiset::swap()` This function is used to exchange the contents of two multisets but the sets must be of same type, although sizes may differ.
	 * lower bound and upper bound for multiset set1
     * `*set1.lower_bound() << endl;`
     * `*set1.upper_bound() << endl;`
