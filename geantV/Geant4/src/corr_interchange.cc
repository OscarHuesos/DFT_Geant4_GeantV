
#include "corr_interchange.hh"
using namespace std;


Correlation_Interchange::Correlation_Interchange(){
}


vector<vector<G4double>>  Correlation_Interchange::lebedev_weights ( G4int siz){

vector<vector<G4double>> leb_matr;
vector<G4double> buffer;

switch (siz) {

case 6:

buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

//cout << "Monday";
break;

case 38:

buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


//10 antes


buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4597);
buffer.push_back(sqrt(1.0-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4597);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4597);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4597);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(0.4597);
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1- 0.4597*0.4597));
buffer.push_back(0.4597);
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

//20

buffer.push_back(sqrt(1 - 0.4597*0.4597));
buffer.push_back(-0.4597);
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(-0.4597);
buffer.push_back(0);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4597);
buffer.push_back(0);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4597);
buffer.push_back(0);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(0.4597);
buffer.push_back(0);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

//25

buffer.push_back(-0.4597);
buffer.push_back(0);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(-0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(0);
buffer.push_back(-0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

//30

buffer.push_back(0);
buffer.push_back(0.4597);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.4597);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.4597);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.4597);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(sqrt(1-0.4597*0.4597));
buffer.push_back(-0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-sqrt(1-0.4597*0.4597));
buffer.push_back(-0.4597);
buffer.push_back(1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();


break;

//////////////////////////////////////////////////////////////////////////////////////
case 86:

buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

//6

buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

//10 ATRAS

buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(-1.0/sqrt(3.0));
buffer.push_back(0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();













// 20 ATRAS

buffer.push_back(0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.3696028);
buffer.push_back(-0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

// 30 ATRAS

buffer.push_back(0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.8525183);
buffer.push_back(0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.8525183);
buffer.push_back(-0.3696028);
buffer.push_back(-0.3696028);
buffer.push_back(0.0111105557);

leb_matr.push_back(buffer);
buffer.clear();

//l2, l2, m2=0.189
//l2=0.694
// 38 ATRAS

buffer.push_back(0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

// 40 atras

buffer.push_back(0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.694354);
buffer.push_back(0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.694354);
buffer.push_back(-0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(-0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

// 50 ATRAS

buffer.push_back(0.694354);
buffer.push_back(0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.694354);
buffer.push_back(0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.694354);
buffer.push_back(-0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1890635);
buffer.push_back(0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1890635);
buffer.push_back(0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1890635);
buffer.push_back(0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

// 60 ATRAS

buffer.push_back(0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.1890635);
buffer.push_back(-0.694354);
buffer.push_back(-0.694354);
buffer.push_back(0.01187650123);

leb_matr.push_back(buffer);
buffer.clear();

// 62 ATRAS
//q1, p1=beta1 (aqui le pusite b1 no se porque en tesis), 0
//p1= 0.927

buffer.push_back(0.374243);
buffer.push_back(0.92733);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.374243);
buffer.push_back(0.92733);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.374243);
buffer.push_back(-0.92733);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.374243);
buffer.push_back(-0.92733);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.92733);
buffer.push_back(0.374243);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.92733);
buffer.push_back(0.374243);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.92733);
buffer.push_back(-0.374243);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.92733);
buffer.push_back(-0.374243);
buffer.push_back(0);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.374243);
buffer.push_back(0);
buffer.push_back(0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.374243);
buffer.push_back(0);
buffer.push_back(0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.374243);
buffer.push_back(0);
buffer.push_back(-0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.374243);
buffer.push_back(0);
buffer.push_back(-0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.92733);
buffer.push_back(0);
buffer.push_back(0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.92733);
buffer.push_back(0);
buffer.push_back(0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.92733);
buffer.push_back(0);
buffer.push_back(-0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.92733);
buffer.push_back(0);
buffer.push_back(-0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.374243);
buffer.push_back(0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.374243);
buffer.push_back(0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

// 80 ATRAS

buffer.push_back(0);
buffer.push_back(0.374243);
buffer.push_back(-0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(0);
buffer.push_back(-0.374243);
buffer.push_back(-0.92733);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.92733);
buffer.push_back(0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.92733);
buffer.push_back(0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.92733);
buffer.push_back(-0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.92733);
buffer.push_back(-0.374243);
buffer.push_back(0.0118123037);

leb_matr.push_back(buffer);
buffer.clear();


break;
case 194:

//a referencia a coeficentes del generador de la pagina:
//https://appliedacousticschalmers.github.io/sound_field_analysis-py/_modules/sound_field_analysis/lebedev.html

buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(1.0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0);
buffer.push_back(-1.0);
buffer.push_back(0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(0);
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

//10 atras

buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

//15
buffer.push_back(-1.0/(sqrt(2) ));
buffer.push_back(1.0/(sqrt(2) ));
buffer.push_back(0);
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0);
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0);
buffer.push_back(0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

//18
//aaa

buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

//20

buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(-1.0/(sqrt(2.0) ));
buffer.push_back(0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

//26
// aab
//segundo valor p1: 0.6713
//q1=p1^2
//lk=sqrt( 1-2 q1*q1)

buffer.push_back(0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

//30
buffer.push_back(0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

//34 atras

buffer.push_back(0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

// 40 atras

buffer.push_back(0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.6712973);
buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.6712973*0.6712973));
buffer.push_back(-0.6712973);
buffer.push_back(-0.6712973);
buffer.push_back(0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

//tal vez tercer valor p1^2: 0.289
//q1=p1^2
//lk=sqrt( 1-2 q1*q1) =

// aab(0.00516, 0.289)
//cof a, ob = 0.5158237711805383e-2
// 50 atras

buffer.push_back(0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

// 60 atras

buffer.push_back(0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.2892466);
buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();









buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

// 70 atras

buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-sqrt(1.0 - 2.0*0.2892466*0.2892466));
buffer.push_back(-0.2892466);
buffer.push_back(-0.2892466);
buffer.push_back(0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

// 74 atras
//aab ( p1^2 = 0.4447  ,conf = 0.005518)


buffer.push_back(0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

// 80 atras

buffer.push_back(0.4447);
buffer.push_back(-0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(-0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(-0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.4447);
buffer.push_back(-0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.4447);
buffer.push_back(-0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

// 90 atras:

buffer.push_back(0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.77749);
buffer.push_back(0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.77749);
buffer.push_back(-0.4447);
buffer.push_back(0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.77749);
buffer.push_back(0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.77749);
buffer.push_back(0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.77749);
buffer.push_back(-0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.77749);
buffer.push_back(-0.4447);
buffer.push_back(-0.4447);
buffer.push_back(0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

// 98 atras
// aab:
//p1^2 = 0.13, a, C_f, lk o mk = 0.4106777028169394e-2)
//no se o como se pongan los coff aqui )

//b = sqrt(1-2aa) =   sqrt(1-2 p1^2 p1^2 ) ? = 0.983

buffer.push_back(0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.98297);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

// 100 atras

buffer.push_back(0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1299335);
buffer.push_back(0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(-0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

// 110

buffer.push_back(0.1299335);
buffer.push_back(0.983);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(0.983);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.1299335);
buffer.push_back(-0.983);
buffer.push_back(-0.1299335);
buffer.push_back(00.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.1299335);
buffer.push_back(-0.983);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.983);
buffer.push_back(0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.983);
buffer.push_back(-0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.983);
buffer.push_back(-0.1299335);
buffer.push_back(0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.983);
buffer.push_back(0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.983);
buffer.push_back(0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

// 120 atras

buffer.push_back(0.983);
buffer.push_back(-0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.983);
buffer.push_back(-0.1299335);
buffer.push_back(-0.1299335);
buffer.push_back(0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

// 122 atras
// ab0:
//(0.5051846064614808e-2, 0.3457702197611283e0)
//p1^2 = 0.3457, a, coff = 0.00505)

//b = sqrt(1-2aa) =   sqrt(1- p1^2 p1^2 ) ?
// = sqrt(1 - 0.3457*0.3457) aprox 0.938

buffer.push_back(0.34577);
buffer.push_back(0.93832);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.34577);
buffer.push_back(0.93832);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.34577);
buffer.push_back(-0.93832);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.34577);
buffer.push_back(-0.93832);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.93832);
buffer.push_back(0.34577);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();


buffer.push_back(-0.93832);
buffer.push_back(0.34577);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.93832);
buffer.push_back(-0.34577);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.93832);
buffer.push_back(-0.34577);
buffer.push_back(0);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

// 130 atras

buffer.push_back(0.34577);
buffer.push_back(0);
buffer.push_back(0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.34577);
buffer.push_back(0);
buffer.push_back(0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.34577);
buffer.push_back(0);
buffer.push_back(-0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.34577);
buffer.push_back(0);
buffer.push_back(-0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.93832);
buffer.push_back(0);
buffer.push_back(0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.93832);
buffer.push_back(0);
buffer.push_back(0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.93832);
buffer.push_back(0);
buffer.push_back(-0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.93832);
buffer.push_back(0);
buffer.push_back(-0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.34577);
buffer.push_back(0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.34577);
buffer.push_back(0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

// 140 atras

buffer.push_back(0);
buffer.push_back(0.34577);
buffer.push_back(-0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.34577);
buffer.push_back(-0.93832);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.93832);
buffer.push_back(0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.93832);
buffer.push_back(0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(0.93832);
buffer.push_back(-0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0);
buffer.push_back(-0.93832);
buffer.push_back(-0.34577);
buffer.push_back(0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

//abc:
//0.5530248916233094e-2, 0.1590417105383530e0, 0.8360360154824589e0)
//p1^2 (1), lk(1) o mk(1) segun corresponda checar = 0.16, a, coff = 0.00553)
//p2^2 (2), lk(2) o mk(2) segun corresponda checar = 0.836
//c = sqrt(1.0 - a*a - b*b) = sqrt(1.0 - 0.16*0.16 - 0.836*0.836 ) aprox 0.525
// 146 atras

buffer.push_back(0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 150 atras

buffer.push_back(0.159042);
buffer.push_back(0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.159042);
buffer.push_back(-0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(-0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.159042);
buffer.push_back(0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 160 atras

buffer.push_back(0.159042);
buffer.push_back(-0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.159042);
buffer.push_back(-0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(-0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(-0.159042);
buffer.push_back(-0.52512);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 170 atras 

buffer.push_back(0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(0.52512);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(-0.52512);
buffer.push_back(0.16);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.836036);
buffer.push_back(-0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.836036);
buffer.push_back(-0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 180

buffer.push_back(0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(-0.159042);
buffer.push_back(0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(-0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(-0.159042);
buffer.push_back(-0.836036);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(-0.836036);
buffer.push_back(0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 190 atras

buffer.push_back(0.52512);
buffer.push_back(0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(0.52512);
buffer.push_back(-0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

buffer.push_back(-0.52512);
buffer.push_back(-0.836036);
buffer.push_back(-0.159042);
buffer.push_back(0.005530249);

leb_matr.push_back(buffer);
buffer.clear();


break;
}

return leb_matr;
}


G4double Correlation_Interchange::evaluation_function( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){


G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);


return N*(pow(x - X_At, L[0]))*(pow(y - Y_At, L[1]))*(pow(z - Z_At, L[2]) )*G4Exp(-alpha*rA2);
}



G4double Correlation_Interchange::evaluation_partial_function_x( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){
 
G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (x - X_At) == 0  ){
x = x + 0.000001; 
}

G4double res_borrar = N*( G4Exp(-alpha*rA2)*(  L[0]*pow(x - X_At, L[0]-1)*(pow(y - Y_At, L[1]))*(pow(z - Z_At, L[2]) ) -
2.0*alpha*(pow(x - X_At , L[0] + 1 ))*(pow(y - Y_At, L[1]))*(pow(z - Z_At, L[2]) )       )  );


if(  isnan( res_borrar ) == true  ) { 
printf("is n l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f N %f rA2 %f alpha %f \n", 
L[0] , L[1] , L[2] , x, y, z , X_At, Y_At,  Z_At ,N , rA2 , alpha  );

}

return res_borrar;

}


G4double Correlation_Interchange::evaluation_partial_function_y( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){


G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}


return N*( G4Exp(-alpha*rA2)*(  L[1]*pow(x - X_At, L[0])*(pow(y - Y_At, L[1] - 1 ))*(pow(z - Z_At, L[2]) ) -
2.0*alpha*(pow(x - X_At , L[0]  ))*(pow(y - Y_At, L[1] + 1 ))*(pow(z - Z_At, L[2]) )       )  );
}


G4double Correlation_Interchange::evaluation_partial_function_z( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){


G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}


return N*( G4Exp(-alpha*rA2)*( L[2]*pow(x - X_At, L[0])*(pow(y - Y_At, L[1]  ))*(pow(z - Z_At, L[2] - 1) ) -
2.0*alpha*(pow(x - X_At , L[0]  ))*(pow(y - Y_At, L[1] ))*(pow(z - Z_At, L[2] + 1) )    )  );
}


G4double Correlation_Interchange::evaluation_double_partial_function_x( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);


if( (x - X_At) == 0  ){
x = x + 0.000001; 
}

G4double ev = N*( G4Exp(-alpha*rA2)*pow(y - Y_At, L[1] )*pow(z - Z_At , L[2])*
( ( L[0]*L[0] - L[0] )*pow(x - X_At, L[0] - 2) - 2*alpha*pow(x - X_At, L[0] )*(2*L[0] + 1) + 4*alpha*alpha*pow( x - X_At , L[0] + 2 )  )  );

return ev;
}


G4double Correlation_Interchange::evaluation_double_partial_function_y( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}

// necesario N
return N*( G4Exp(-alpha*rA2)*pow(x - X_At, L[0] )*pow(z - Z_At , L[2])*
( ( L[1]*L[1] - L[1] )*pow(y - Y_At, L[1] - 2) 
-2*alpha*pow(y - Y_At, L[1] )*(2*L[1] + 1) + 4*alpha*alpha*pow( y - Y_At , L[1] + 2 )  )  );
}



G4double Correlation_Interchange::evaluation_double_partial_function_z( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow(x - X_At, L[0] )*pow(y - Y_At , L[1])*
( ( L[2]*L[2] - L[2] )*pow( z - Z_At, L[2] - 2) 
-2*alpha*pow(z - Z_At, L[2] )*(2*L[2] + 1) + 4*alpha*alpha*pow( z - Z_At , L[2] + 2 )  )  );
}



G4double Correlation_Interchange::evaluation_partial_function_dydx( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (x - X_At) == 0  ){
x = x + 0.000001; 
}


if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow(z - Z_At , L[2])*( L[0]*pow( x - X_At , L[0] - 1)
-2*alpha*pow(x - X_At, L[0] + 1) )*( L[1]*pow( y - Y_At , L[1] - 1) -2*alpha*pow(y - Y_At, L[1] + 1) ) );
}


G4double Correlation_Interchange::evaluation_partial_function_dzdx( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}

if( (x - X_At) == 0  ){
x = x + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow(y - Y_At , L[1])*( L[0]*pow( x - X_At , L[0] - 1)
-2*alpha*pow(x - X_At, L[0] + 1) )*( L[2]*pow( z - Z_At , L[2] - 1) -2*alpha*pow(z - Z_At, L[2] + 1) ) );
}



G4double Correlation_Interchange::evaluation_partial_function_dxdy( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);


if( (x - X_At) == 0  ){
x = x + 0.000001; 
}

if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*(  L[1]*pow( y - Y_At, L[1] -1 ) -
2*alpha*pow( y - Y_At, L[1] + 1) )*pow(z - Z_At, L[2])*(  L[0]*pow( x - X_At , L[0] -1) - 2*alpha*pow( x - X_At , L[0] + 1) ) );
}


G4double Correlation_Interchange::evaluation_partial_function_dxdz( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (x - X_At) == 0  ){
x = x + 0.000001; 
}

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow( y - Y_At, L[1] )*( L[2]*pow( z - Z_At , L[2] -1 ) -2*alpha*pow( z - Z_At , L[2] + 1)  )*
(  L[0]*pow(x - X_At , L[0] - 1)  -2*alpha*pow( x - X_At , L[0] + 1)   ) );
}


G4double Correlation_Interchange::evaluation_partial_function_dzdy( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);

if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow( x - X_At , L[0])*( L[1]*pow( y - Y_At , L[1] - 1) - 
2*alpha*pow(y - Y_At , L[1] + 1) )*( L[2]*pow( z - Z_At , L[2] - 1) - 2*alpha*pow( z - Z_At , L[2] + 1)  )  );
}



G4double Correlation_Interchange::evaluation_partial_function_dydz( G4double X_At,
G4double Y_At, G4double Z_At, G4double alpha, G4double N,
vector<G4int> L, G4double x, G4double y, G4double z ){

G4double rA2 = (x - X_At)*(x - X_At) + (y - Y_At)*(y - Y_At) +
(z - Z_At)*(z - Z_At);


if( (y - Y_At) == 0  ){
y = y + 0.000001; 
}

if( (z - Z_At) == 0  ){
z = z + 0.000001; 
}

return N*( G4Exp(-alpha*rA2)*pow( x - X_At , L[0])*( L[2]*pow( z - Z_At , L[2] -1 ) -
2*alpha*pow( z - Z_At , L[2] + 1) )*( L[1]*pow( y - Y_At, L[1] - 1) -2*alpha*pow( y - Y_At , L[1] + 1)  ) );
}




vector<G4double> Correlation_Interchange::becke_eval ( Atom& atomo_puntos ,  Molecule& Mol, G4int& cont , G4int id , G4int functional ){

vector<Residue *>::iterator res_i;
vector<Atom *>::iterator atr_i;

vector<Residue *>::iterator res_j;
vector<Atom *>::iterator atr_j;

vector<G4double> VBpoint;
vector<G4double> pesos;

vector<Residue *>::iterator res_k;
vector<Atom *>::iterator atr_k;

G4int size_puntos = atomo_puntos.atomic_grid.size();
G4int vpuntos;
G4double Var = 1.0;
G4double x,y,z;
G4double fk, mu, ra, rb, R, uab, acc, a, wa;
G4double  bb, bx, by, bz;
G4double ev = 0, evx = 0, evy = 0, evz = 0;
G4int dd = 0;

G4double dxdx, dydy, dzdz;
G4double dxdy, dxdz;
G4double dydx, dydz;
G4double dzdx, dzdy;

G4double sdxdx = 0, sdydy = 0, sdzdz = 0;
G4double sdxdy = 0, sdxdz = 0;
G4double sdydx = 0, sdydz = 0;
G4double sdzdx = 0, sdzdy = 0;

for(int i = 0; i < size_puntos ; i++) {

G4int size_malla = atomo_puntos.atomic_grid[i].leb.size();

for(int j = 0; j < size_malla ; j++) {

x = atomo_puntos.atomic_grid[i].leb[j][0];
y = atomo_puntos.atomic_grid[i].leb[j][1];
z = atomo_puntos.atomic_grid[i].leb[j][2];

G4int atm_cont = 0;

for(res_i = Mol.Residuos_cadena.begin();
    res_i != Mol.Residuos_cadena.end() ; ++res_i ){

for(atr_i = (*res_i)->Lista_atoms.begin();
     atr_i != (*res_i)->Lista_atoms.end(); ++ atr_i ){

ra =  sqrt(  ((*atr_i)->fX - x)*((*atr_i)->fX - x) +
((*atr_i)->fY - y)*((*atr_i)->fY - y) + ((*atr_i)->fZ - z)*
((*atr_i)->fZ - z) );


dd++;
atm_cont++;

for(res_j = Mol.Residuos_cadena.begin();
    res_j != Mol.Residuos_cadena.end() ; ++res_j ){

 for(atr_j = (*res_j)->Lista_atoms.begin();
     atr_j != (*res_j)->Lista_atoms.end(); ++ atr_j ){

   if (  ( (*atr_i)->fX  == (*atr_j)->fX ) &&
         ( (*atr_i)->fY  == (*atr_j)->fY  ) &&
         ( (*atr_i)->fZ  == (*atr_j)->fZ)  ){


}else{

rb = sqrt(  ((*atr_j)->fX - x)*((*atr_j)->fX - x) +
((*atr_j)->fY - y)*((*atr_j)->fY - y) + ((*atr_j)->fZ - z)*
((*atr_j)->fZ - z) );

R = sqrt(  ( (*atr_i)->fX - (*atr_j)->fX)*( (*atr_i)->fX - (*atr_j)->fX) +
( (*atr_i)->fY - (*atr_j)->fY)*( (*atr_i)->fY - (*atr_j)->fY) +
( (*atr_i)->fZ - (*atr_j)->fZ)*( (*atr_i)->fZ - (*atr_j)->fZ) );


mu = (ra-rb)/(R);

//if heteronuclear
if(  (*atr_i)->fZat  !=   (*atr_j)->fZat  ){

uab = ( (*atr_i)->Bragg_radius -
(*atr_j)->Bragg_radius )/((*atr_i)->Bragg_radius +
(*atr_j)->Bragg_radius );

a = uab/((uab*uab)-1);

//modfy bt 0.45
if(a  > 0.5  ){
a =  0.5;
}

if(a  < -0.5  ){
a = -0.5;
}


fk = mu + a*(1-mu*mu);

// if homonuclear
}else{
fk = mu;
}


for(int ii = 0; ii < 3 ; ii++) {
fk = fk*(3.0 - fk*fk)/2.0;
}


Var = Var*(1 - fk )/2.0;

}

}
}


VBpoint.push_back(Var);
Var = 1.0;
}
} //closing ij atomic

acc = 0;
vpuntos = VBpoint.size();

for(int p = 0; p < vpuntos ; p++) {
acc = acc + VBpoint[p];
}

wa = VBpoint[id]/acc;
pesos.push_back(atomo_puntos.atomic_grid[i].leb[j][3]*wa);
VBpoint.clear();

G4int cont_int = 0;

for(res_k = Mol.Residuos_cadena.begin();
    res_k != Mol.Residuos_cadena.end() ; ++res_k){

 for(atr_k = (*res_k)->Lista_atoms.begin();
     atr_k != (*res_k)->Lista_atoms.end(); ++ atr_k){
 
  for((*atr_k)->shell1 = (*atr_k)->orbitals.begin();
      (*atr_k)->shell1 != (*atr_k)->orbitals.end(); ++ (*atr_k)->shell1){

for (int k = 0; k < 3; k++) {


bb = evaluation_function( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);


// if gga, partial derivates respect to thr position:

ev = ev + (  (*(*atr_k)->shell1).NAC[k].C )*bb;


if(functional != 0 ){

bx = evaluation_partial_function_x ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);


if(  isnan(bx)   ) { 
printf("bx is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

evx = evx + (  (*(*atr_k)->shell1).NAC[k].C )*bx;

by = evaluation_partial_function_y ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

evy = evy + (  (*(*atr_k)->shell1).NAC[k].C )*by;

bz = evaluation_partial_function_z ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

evz = evz + (  (*(*atr_k)->shell1).NAC[k].C )*bz;

// double derivatives

dxdx = evaluation_double_partial_function_x ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);


if(  isnan(dxdx)   ) { 
printf("dxdx is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdxdx = sdxdx + (  (*(*atr_k)->shell1).NAC[k].C )*dxdx;

dydy = evaluation_double_partial_function_y ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dydy)   ) { 
printf("dydy is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdydy = sdydy + ( (*(*atr_k)->shell1).NAC[k].C )*dydy;

dzdz = evaluation_double_partial_function_z ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dzdz)   ) { 
printf("dzdz is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdzdz = sdzdz + ( (*(*atr_k)->shell1).NAC[k].C )*dzdz;

dxdy =  evaluation_partial_function_dxdy ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dxdy)   ) { 
printf("dxdy is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}


sdxdy = sdxdy + ( (*(*atr_k)->shell1).NAC[k].C )*dxdy; 

dxdz = evaluation_partial_function_dxdz ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dxdz)   ) { 
printf("dxdz is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}


sdxdz = sdxdz + ( (*(*atr_k)->shell1).NAC[k].C )*dxdz;

dydx = evaluation_partial_function_dydx ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dydx)   ) { 
printf("dydx is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdydx = sdydx + ( (*(*atr_k)->shell1).NAC[k].C )*dydx;

dydz = evaluation_partial_function_dydz ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dydz)   ) { 
printf("dydz is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdydz = sdydz + ( (*(*atr_k)->shell1).NAC[k].C )*dydz;

dzdx = evaluation_partial_function_dzdx ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dzdx)   ) { 
printf("dzdx is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdzdx = sdzdx + ( (*(*atr_k)->shell1).NAC[k].C )*dzdx;

dzdy = evaluation_partial_function_dzdy ( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
(*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
(*(*atr_k)->shell1).l, x ,y, z);

if(  isnan(dzdy)   ) { 
printf("dzdy is nan with l1 %d l2 %d l3 %d x %f y %f z %f Ax %f Ay %f Az %f \n",  (*(*atr_k)->shell1).l[0], (*(*atr_k)->shell1).l[1],
 (*(*atr_k)->shell1).l[2],  x , y, z, (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ  );
}

sdzdy = sdzdy + ( (*(*atr_k)->shell1).NAC[k].C )*dzdy;

}

}


(*(*atr_k)->shell1).Func.push_back(ev);

if(functional != 0 ){


(*(*atr_k)->shell1).x_grad.push_back(evx);
(*(*atr_k)->shell1).y_grad.push_back(evy);
(*(*atr_k)->shell1).z_grad.push_back(evz);

(*(*atr_k)->shell1).dxx.push_back(sdxdx);
(*(*atr_k)->shell1).dyy.push_back(sdydy);
(*(*atr_k)->shell1).dzz.push_back(sdzdz);

(*(*atr_k)->shell1).dxy.push_back(sdxdy);
(*(*atr_k)->shell1).dxz.push_back(sdxdz);

(*(*atr_k)->shell1).dyx.push_back(sdydx);
(*(*atr_k)->shell1).dyz.push_back(sdydz);

(*(*atr_k)->shell1).dzx.push_back(sdzdx);
(*(*atr_k)->shell1).dzy.push_back(sdzdy);

}

ev = 0;
evx=0;
evy=0;
evz=0;

sdxdx=0;
sdydy=0;
sdzdz=0;

sdxdy = 0;
sdxdz = 0;

sdydx = 0;
sdydz = 0;

sdzdx = 0;
sdzdy = 0;

cont_int++;
}
}
}


cont++;
}// j
}// i

return pesos;
}


void Correlation_Interchange::Euler_Mac_Lebedev_SG1 ( Atom& atomo){

vector<vector<G4double>>::iterator Col;
vector<G4double>::iterator Row;

G4int Nr = 50;
G4double wir =0;
G4double ri = 0;
G4double bf = 0;
G4int No = 0;
G4double RE = 0;
G4int cont = 0;
/////////////////////////////////////////////////////

//sg-0

//sg-1

for (int i = 1; i < Nr + 1 ; i++) {

EML points;
bf = Nr+1.0-i;

wir = 2*atomo.Bohr_radius*atomo.Bohr_radius*atomo.Bohr_radius*
( (Nr+1.0)*pow(i,5))/(pow(bf,7) );

ri = atomo.Bohr_radius*(i*i/( bf*bf ));

if( ( atomo.Bohr_radius == 1.0  ) || ( atomo.Bohr_radius == 4.0/7.0 )   ){

RE= (0.25)*atomo.Bohr_radius;

if( ri < RE){
No = 6;
}else{
RE= (0.5)*atomo.Bohr_radius;
if( ri < RE){
No = 38;
}else{
RE = atomo.Bohr_radius;
if (ri < RE){
No = 86;
}else{
RE = (9.0/2.0)*atomo.Bohr_radius;
if( ri < RE){
No = 194;
}else{
// default

No = 86;
}
}
}
}

}else{

if( ( atomo.Bohr_radius ==  9.0/7.0) || ( atomo.Bohr_radius ==  1.0563)  ||  (  atomo.Bohr_radius ==  0.9)   ){

RE= 0.1667*atomo.Bohr_radius;

if( ri < RE){
No = 6;
}else{
RE= (0.5)*atomo.Bohr_radius;
if( ri < RE){
No = 38;
}else{
RE = 0.9*atomo.Bohr_radius;
if (ri < RE){
No = 86;
}else{
RE = (7.0/2.0)*atomo.Bohr_radius;
if( ri < RE){
No = 194;
}else{

No = 86;
}
}
}
}

}else{

if( ( atomo.Bohr_radius == 1.7862) || ( atomo.Bohr_radius ==  1.6489)  ){

RE = 0.1*atomo.Bohr_radius;

if( ri < RE){
No = 6;
}else{
RE= 0.4*atomo.Bohr_radius;
if( ri < RE){
No = 38;
}else{
RE = 0.8*atomo.Bohr_radius;
if (ri < RE){
No = 86;
}else{
RE = (5.0/2.0)*atomo.Bohr_radius;
if( ri < RE){
No = 194;
}else{

//////////////////////////////////////////////////////////////////
No = 86;
}
}
}
}

}
}
}



points.Nomega = No;
points.leb  = lebedev_weights (No);
G4int siz =  points.leb.size();

for(int i = 0; i <  siz ; i++) {

points.leb[i][0] =  points.leb[i][0]*ri + atomo.fX;
points.leb[i][1] =  points.leb[i][1]*ri + atomo.fY;
points.leb[i][2] =  points.leb[i][2]*ri + atomo.fZ;
points.leb[i][3] =  points.leb[i][3]*wir*(4.0)*pi;

cont++;
}

atomo.atomic_grid.push_back(points);

}

return;
}



vector<G4double>  Correlation_Interchange::Grid( Molecule& Mol , G4int funcional ){

vector<vector<G4double>>::iterator col_puntos;


vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;

vector<vector<G4double>>::iterator Col;
vector<G4double>  W;
vector<G4double> buf_Sum;
vector<EML >::iterator eml;

G4int sis = 0;
G4int cont = 0;

for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

  for(atr_main = (*res_main)->Lista_atoms.begin();
      atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

Euler_Mac_Lebedev_SG1 ( (*(*atr_main))  );

}
}

Mol.No_points = 0;
for(res_main = Mol.Residuos_cadena.begin();
    res_main != Mol.Residuos_cadena.end() ; ++res_main){

 for(atr_main = (*res_main)->Lista_atoms.begin();
    atr_main != (*res_main)->Lista_atoms.end(); ++ atr_main){

sis = 0;
buf_Sum = becke_eval (  (*(*atr_main)) ,  Mol , sis, cont, funcional );
Mol.No_points = Mol.No_points + sis;

W.insert(W.end(), buf_Sum.begin(), buf_Sum.end());
buf_Sum.clear();


cont++;
}
}

return W;
}


