#include <iostream>

#include<array3dbox.h>

using namespace std;
int main(){
	Array3dbox<int> a(2);
	
	a.push_element(0.1, 0.0, 0.0, 33);
	
	a.push_element(-0.1, 0.0, 0.0, 23);
	
	for(auto p=a.begin(); p!=a.end(); ++p){
		cout << p->first << "; " ;
		vector<int> &v=p->second;
		for(auto i=v.begin(); i!=v.end(); ++i)
			cout << *i << ' ';
			
		cout << '\n';
	}
	return 0;
}
