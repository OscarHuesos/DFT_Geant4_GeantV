

#include "corr_interchange.hh"
using namespace std;

Correlation_Interchange::Correlation_Interchange(){
}


vector<vector<Double_v>>  Correlation_Interchange::lebedev_weightsV ( G4int siz){

vector<vector<Double_v>> lebedev_total;
vector<Double_v> lebedev_buffer;

vector<vector<G4double>> leb_matr;
vector<G4double> buffer;

Double_v buvx;
Double_v buvy;
Double_v buvz;
Double_v buvw;

switch (siz) {

case 6:

Set(buvx, 0, 1.0);
Set(buvy, 0, 0);
Set(buvz, 0, 0);
Set(buvw, 0, 1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0);
Set(buvy, 1, 0);
Set(buvz, 1, 0);
Set(buvw, 1, 1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 2, 0);
Set(buvy, 2, 1.0);
Set(buvz, 2, 0);
Set(buvw, 2, 1.0/6.0);


leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 3, 0);
Set(buvy, 3, -1.0);
Set(buvz, 3, 0);
Set(buvw, 3, 1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();


lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 0);
Set(buvz, 0, 1.0);
Set(buvw, 0, 1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, 0);
Set(buvz, 1, -1.0);
Set(buvw, 1, 1.0/6.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0);
Set(buvz, 2, 0);
Set(buvw, 2, 0);

Set(buvx, 3, 0);
Set(buvy, 3, 0);
Set(buvz, 3, 0);
Set(buvw, 3, 0);


lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

break;

case 38:

Set(buvx, 0, 1.0);
Set(buvy, 0, 0);
Set(buvz, 0, 0);
Set(buvw, 0, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0);
Set(buvy, 1, 0);
Set(buvz, 1, 0);
Set(buvw, 1, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 1.0);
Set(buvz, 2, 0);
Set(buvw, 2, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -1.0);
Set(buvz, 3, 0);
Set(buvw, 3, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();


lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 0);
Set(buvz, 0, 1.0);
Set(buvw, 0, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 1, 0);
Set(buvy, 1, 0);
Set(buvz, 1, -1.0);
Set(buvw, 1, 1.0/105.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 1.0/sqrt(3.0));
Set(buvy, 2, 1.0/sqrt(3.0));
Set(buvz, 2, 1.0/sqrt(3.0));
Set(buvw, 2, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/sqrt(3.0));
Set(buvy, 3, 1.0/sqrt(3.0));
Set(buvz, 3, 1.0/sqrt(3.0));
Set(buvw, 3, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/sqrt(3.0));
Set(buvy, 0, -1.0/sqrt(3.0));
Set(buvz, 0, 1.0/sqrt(3.0));
Set(buvw, 0, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/sqrt(3.0));
Set(buvy, 1, -1.0/sqrt(3.0));
Set(buvz, 1, 1.0/sqrt(3.0));
Set(buvw, 1, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 1.0/sqrt(3.0));
Set(buvy, 2, 1.0/sqrt(3.0));
Set(buvz, 2, -1.0/sqrt(3.0));
Set(buvw, 2, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/sqrt(3.0));
Set(buvy, 3, 1.0/sqrt(3.0));
Set(buvz, 3, -1.0/sqrt(3.0));
Set(buvw, 3, 9.0/280.0);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/sqrt(3.0));
Set(buvy, 0, -1.0/sqrt(3.0));
Set(buvz, 0, -1.0/sqrt(3.0));
Set(buvw, 0, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/sqrt(3.0));
Set(buvy, 1, -1.0/sqrt(3.0));
Set(buvz, 1, -1.0/sqrt(3.0));
Set(buvw, 1, 9.0/280.0);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 2, 0.46);
Set(buvy, 2, sqrt(1-0.46*0.46) );
Set(buvz, 2, 0);
Set(buvw, 2, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 3, -0.46);
Set(buvy, 3, sqrt(1-0.46*0.46) );
Set(buvz, 3, 0);
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.46);
Set(buvy, 0, -sqrt(1-0.46*0.46) );
Set(buvz, 0, 0);
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.46);
Set(buvy, 1, -sqrt(1-0.46*0.46) );
Set(buvz, 1, 0);
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, sqrt(1-0.46*0.46) );
Set(buvy, 2, 0.46 );
Set(buvz, 2, 0);
Set(buvw, 2, 1.0/35.0);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1-0.46*0.46) );
Set(buvy, 3, 0.46 );
Set(buvz, 3, 0);
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

//20

Set(buvx, 0, sqrt(1-0.46*0.46) );
Set(buvy, 0, -0.46 );
Set(buvz, 0, 0);
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1-0.46*0.46) );
Set(buvy, 1, -0.46 );
Set(buvz, 1, 0);
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.46); 
Set(buvy, 2, 0);
Set(buvz, 2, sqrt(1-0.46*0.46) );
Set(buvw, 2, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.46); 
Set(buvy, 3, 0);
Set(buvz, 3, sqrt(1-0.46*0.46) );
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();


lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.46); 
Set(buvy, 0, 0);
Set(buvz, 0, -sqrt(1-0.46*0.46) );
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.46); 
Set(buvy, 1, 0);
Set(buvz, 1, -sqrt(1-0.46*0.46) );
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, sqrt(1-0.46*0.46)); 
Set(buvy, 2, 0);
Set(buvz, 2, 0.46 );
Set(buvw, 2, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1-0.46*0.46)); 
Set(buvy, 3, 0);
Set(buvz, 3, 0.46 );
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, sqrt(1-0.46*0.46)); 
Set(buvy, 0, 0);
Set(buvz, 0, -0.46 );
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1-0.46*0.46)); 
Set(buvy, 1, 0);
Set(buvz, 1, -0.46 );
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

//30

Set(buvx, 2, 0  ); 
Set(buvy, 2, 0.46);
Set(buvz, 2, sqrt(1-0.46*0.46) );
Set(buvw, 2, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0  ); 
Set(buvy, 3, -0.46);
Set(buvz, 3, sqrt(1-0.46*0.46) );
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0  ); 
Set(buvy, 0, 0.46);
Set(buvz, 0, -sqrt(1-0.46*0.46) );
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0  ); 
Set(buvy, 1, -0.46);
Set(buvz, 1, -sqrt(1-0.46*0.46) );
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0 ); 
Set(buvy, 2, sqrt(1-0.46*0.46));
Set(buvz, 2, 0.46 );
Set(buvw, 2, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0 ); 
Set(buvy, 3, -sqrt(1-0.46*0.46));
Set(buvz, 3, 0.46 );
Set(buvw, 3, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0 ); 
Set(buvy, 0, sqrt(1-0.46*0.46));
Set(buvz, 0, -0.46 );
Set(buvw, 0, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0 ); 
Set(buvy, 1, -sqrt(1-0.46*0.46));
Set(buvz, 1, -0.46 );
Set(buvw, 1, 1.0/35.0);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0);
Set(buvz, 2, 0);
Set(buvw, 2, 0);

Set(buvx, 3, 0);
Set(buvy, 3, 0);
Set(buvz, 3, 0);
Set(buvw, 3, 0);


lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

break;

//////////////////////////////////////////////////////////////////////////////////////
case 86:


Set(buvx, 0, 1.0);
Set(buvy, 0, 0);
Set(buvz, 0, 0);
Set(buvw, 0, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0);
Set(buvy, 1, 0);
Set(buvz, 1, 0);
Set(buvw, 1, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 1.0);
Set(buvz, 2, 0);
Set(buvw, 2, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -1.0);
Set(buvz, 3, 0);
Set(buvw, 3, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 0);
Set(buvz, 0, 1.0);
Set(buvw, 0, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, 0);
Set(buvz, 1, -1.0);
Set(buvw, 1, 0.0115440115);

leb_matr.push_back(buffer);
buffer.clear();

//6

Set(buvx, 2, 1.0/sqrt(3.0));
Set(buvy, 2, 1.0/sqrt(3.0));
Set(buvz, 2, 1.0/sqrt(3.0));
Set(buvw, 2, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/sqrt(3.0));
Set(buvy, 3, 1.0/sqrt(3.0));
Set(buvz, 3, 1.0/sqrt(3.0));
Set(buvw, 3, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/sqrt(3.0));
Set(buvy, 0, -1.0/sqrt(3.0));
Set(buvz, 0, 1.0/sqrt(3.0));
Set(buvw, 0, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/sqrt(3.0));
Set(buvy, 1, -1.0/sqrt(3.0));
Set(buvz, 1, 1.0/sqrt(3.0));
Set(buvw, 1, 0.011943909);


leb_matr.push_back(buffer);
buffer.clear();

//10 ATRAS

Set(buvx, 2, 1.0/sqrt(3.0));
Set(buvy, 2, 1.0/sqrt(3.0));
Set(buvz, 2, -1.0/sqrt(3.0));
Set(buvw, 2, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/sqrt(3.0));
Set(buvy, 3, 1.0/sqrt(3.0));
Set(buvz, 3, -1.0/sqrt(3.0));
Set(buvw, 3, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, 1.0/sqrt(3.0));
Set(buvy, 0, -1.0/sqrt(3.0));
Set(buvz, 0, -1.0/sqrt(3.0));
Set(buvw, 0, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/sqrt(3.0));
Set(buvy, 1, -1.0/sqrt(3.0));
Set(buvz, 1, -1.0/sqrt(3.0));
Set(buvw, 1, 0.011943909);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.3696028);
Set(buvy, 2, 0.3696028);
Set(buvz, 2, 0.8525183);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.3696028);
Set(buvy, 3, 0.3696028);
Set(buvz, 3, 0.8525183);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.3696028);
Set(buvy, 0, -0.3696028);
Set(buvz, 0, 0.8525183);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.3696028);
Set(buvy, 1, -0.3696028);
Set(buvz, 1, 0.8525183);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.3696028);
Set(buvy, 2, 0.3696028);
Set(buvz, 2, -0.8525183);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.3696028);
Set(buvy, 3, 0.3696028);
Set(buvz, 3, -0.8525183);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 20 ATRAS

Set(buvx, 0, 0.3696028);
Set(buvy, 0, -0.3696028);
Set(buvz, 0, -0.8525183);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.3696028);
Set(buvy, 1, -0.3696028);
Set(buvz, 1, -0.8525183);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 2, 0.3696028);
Set(buvy, 2, 0.8525183);
Set(buvz, 2, 0.3696028);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 3, -0.3696028);
Set(buvy, 3, 0.8525183);
Set(buvz, 3, 0.3696028);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.3696028);
Set(buvy, 0, -0.8525183);
Set(buvz, 0, 0.3696028);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 1, -0.3696028);
Set(buvy, 1, -0.8525183);
Set(buvz, 1, 0.3696028);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 2, 0.3696028);
Set(buvy, 2, 0.8525183);
Set(buvz, 2, -0.3696028);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.3696028);
Set(buvy, 3, 0.8525183);
Set(buvz, 3, -0.3696028);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.3696028);
Set(buvy, 0, -0.8525183);
Set(buvz, 0, -0.3696028);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.3696028);
Set(buvy, 1, -0.8525183);
Set(buvz, 1, -0.3696028);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

// 30 ATRAS

Set(buvx, 2, 0.8525183);
Set(buvy, 2, 0.3696028);
Set(buvz, 2, 0.3696028);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.8525183);
Set(buvy, 3, 0.3696028);
Set(buvz, 3, 0.3696028);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.8525183);
Set(buvy, 0, -0.3696028);
Set(buvz, 0, 0.3696028);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.8525183);
Set(buvy, 1, -0.3696028);
Set(buvz, 1, 0.3696028);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.8525183);
Set(buvy, 2, 0.3696028);
Set(buvz, 2, -0.3696028);
Set(buvw, 2, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.8525183);
Set(buvy, 3, 0.3696028);
Set(buvz, 3, -0.3696028);
Set(buvw, 3, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.8525183);
Set(buvy, 0, -0.3696028);
Set(buvz, 0, -0.3696028);
Set(buvw, 0, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.8525183);
Set(buvy, 1, -0.3696028);
Set(buvz, 1, -0.3696028);
Set(buvw, 1, 0.0111105557 );

leb_matr.push_back(buffer);
buffer.clear();

//l2, l2, m2=0.189
//l2=0.694
// 38 ATRAS

Set(buvx, 2, 0.694354);
Set(buvy, 2, 0.694354);
Set(buvz, 2, 0.1890635);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.694354);
Set(buvy, 3, 0.694354);
Set(buvz, 3, 0.1890635);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 40 atras

Set(buvx, 0, 0.694354);
Set(buvy, 0, -0.694354);
Set(buvz, 0, 0.1890635);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.694354);
Set(buvy, 1, -0.694354);
Set(buvz, 1, 0.1890635);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.694354);
Set(buvy, 2, 0.694354);
Set(buvz, 2, -0.1890635);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.694354);
Set(buvy, 3, 0.694354);
Set(buvz, 3, -0.1890635);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, 0.694354);
Set(buvy, 0, -0.694354);
Set(buvz, 0, -0.1890635);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.694354);
Set(buvy, 1, -0.694354);
Set(buvz, 1, -0.1890635);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.694354);
Set(buvy, 2, 0.1890635);
Set(buvz, 2, 0.694354);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.694354);
Set(buvy, 3, 0.1890635);
Set(buvz, 3, 0.694354);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.694354);
Set(buvy, 0, -0.1890635);
Set(buvz, 0, 0.694354);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.694354);
Set(buvy, 1, -0.1890635);
Set(buvz, 1, 0.694354);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

// 50 ATRAS

Set(buvx, 2, 0.694354);
Set(buvy, 2, 0.1890635);
Set(buvz, 2, -0.694354);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.694354);
Set(buvy, 3, 0.1890635);
Set(buvz, 3, -0.694354);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.694354);
Set(buvy, 0, -0.1890635);
Set(buvz, 0, -0.694354);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.694354);
Set(buvy, 1, -0.1890635);
Set(buvz, 1, -0.694354);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.1890635);
Set(buvy, 2, 0.694354);
Set(buvz, 2, 0.694354);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.1890635);
Set(buvy, 3, 0.694354);
Set(buvz, 3, 0.694354);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.1890635);
Set(buvy, 0, -0.694354);
Set(buvz, 0, 0.694354);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.1890635);
Set(buvy, 1, -0.694354);
Set(buvz, 1, 0.694354);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.1890635);
Set(buvy, 2, 0.694354);
Set(buvz, 2, -0.694354);
Set(buvw, 2, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.1890635);
Set(buvy, 3, 0.694354);
Set(buvz, 3, -0.694354);
Set(buvw, 3, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 60 ATRAS

Set(buvx, 0, 0.1890635);
Set(buvy, 0, -0.694354);
Set(buvz, 0, -0.694354);
Set(buvw, 0, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.1890635);
Set(buvy, 1, -0.694354);
Set(buvz, 1, -0.694354);
Set(buvw, 1, 0.01187650123 );

leb_matr.push_back(buffer);
buffer.clear();

// 62 ATRAS
//p1= 0.927

Set(buvx, 2, 0.374243);
Set(buvy, 2, 0.92733);
Set(buvz, 2, 0);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.374243);
Set(buvy, 3, 0.92733);
Set(buvz, 3, 0);
Set(buvw, 3, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.374243);
Set(buvy, 0, -0.92733);
Set(buvz, 0, 0);
Set(buvw, 0, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.374243);
Set(buvy, 1, -0.92733);
Set(buvz, 1, 0);
Set(buvw, 1, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.92733);
Set(buvy, 2, 0.374243);
Set(buvz, 2, 0);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.92733);
Set(buvy, 3, 0.374243);
Set(buvz, 3, 0);
Set(buvw, 3, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.92733);
Set(buvy, 0, -0.374243);
Set(buvz, 0, 0);
Set(buvw, 0, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.92733);
Set(buvy, 1, -0.374243);
Set(buvz, 1, 0);
Set(buvw, 1, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.374243);
Set(buvy, 2, 0);
Set(buvz, 2, 0.92733);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.374243);
Set(buvy, 3, 0);
Set(buvz, 3, 0.92733);
Set(buvw, 3, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.374243);
Set(buvy, 0, 0);
Set(buvz, 0, -0.92733);
Set(buvw, 0, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.374243);
Set(buvy, 1, 0);
Set(buvz, 1, -0.92733);
Set(buvw, 1, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.92733);
Set(buvy, 2, 0);
Set(buvz, 2, 0.374243);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.92733);
Set(buvy, 3, 0);
Set(buvz, 3, 0.374243);
Set(buvw, 3, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.92733);
Set(buvy, 0, 0);
Set(buvz, 0, -0.374243);
Set(buvw, 0, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.92733);
Set(buvy, 1, 0);
Set(buvz, 1, -0.374243);
Set(buvw, 1, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0.374243);
Set(buvz, 2, 0.92733);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -0.374243);
Set(buvz, 3, 0.92733);
Set(buvw, 3, 0.0118123037 );


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 80 ATRAS

Set(buvx, 0, 0);
Set(buvy, 0, 0.374243);
Set(buvz, 0, -0.92733);
Set(buvw, 0, 0.0118123037 );


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, -0.374243);
Set(buvz, 1, -0.92733);
Set(buvw, 1, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0.92733);
Set(buvz, 2, 0.374243);
Set(buvw, 2, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -0.92733);
Set(buvz, 3, 0.374243);
Set(buvw, 3, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, 0);
Set(buvy, 0, 0.92733);
Set(buvz, 0, -0.374243);
Set(buvw, 0, 0.0118123037 );

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, -0.92733);
Set(buvz, 1, -0.374243);
Set(buvw, 1, 0.0118123037 );


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0);
Set(buvz, 2, 0);
Set(buvw, 2, 0 );

Set(buvx, 3, 0);
Set(buvy, 3, 0);
Set(buvz, 3, 0);
Set(buvw, 3, 0);

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


break;
case 194:

Set(buvx, 0, 1.0);
Set(buvy, 0, 0);
Set(buvz, 0, 0);
Set(buvw, 0, 0.00178234044);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0);
Set(buvy, 1, 0);
Set(buvz, 1, 0);
Set(buvw, 1, 0.00178234044);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 1.0);
Set(buvz, 2, 0);
Set(buvw, 2, 0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -1.0);
Set(buvz, 3, 0);
Set(buvw, 3, 0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 0);
Set(buvz, 0, 1.0);
Set(buvw, 0, 0.00178234044);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, 0);
Set(buvz, 1, -1.0);
Set(buvw, 1, 0.00178234044);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 1.0/(sqrt(2) ));
Set(buvz, 2, 1.0/(sqrt(2) ));
Set(buvw, 2, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -1.0/(sqrt(2) ));
Set(buvz, 3, 1.0/(sqrt(2) ));
Set(buvw, 3, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 1.0/(sqrt(2) ));
Set(buvz, 0, -1.0/(sqrt(2) ));
Set(buvw, 0, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, -1.0/(sqrt(2) ));
Set(buvz, 1, -1.0/(sqrt(2) ));
Set(buvw, 1, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

//10 a

Set(buvx, 2, 1.0/(sqrt(2) ));
Set(buvy, 2, 0 );
Set(buvz, 2, 1.0/(sqrt(2) ));
Set(buvw, 2, 0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/(sqrt(2) ));
Set(buvy, 3, 0 );
Set(buvz, 3, 1.0/(sqrt(2) ));
Set(buvw, 3, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/(sqrt(2) ));
Set(buvy, 0, 0 );
Set(buvz, 0, -1.0/(sqrt(2) ));
Set(buvw, 0, 0.005716906);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/(sqrt(2) ));
Set(buvy, 1, 0 );
Set(buvz, 1, -1.0/(sqrt(2) ));
Set(buvw, 1, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 1.0/(sqrt(2) ));
Set(buvy, 2, 1.0/(sqrt(2) ));
Set(buvz, 2, 0 );
Set(buvw, 2, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

//15

Set(buvx, 3, -1.0/(sqrt(2) ));
Set(buvy, 3, 1.0/(sqrt(2) ));
Set(buvz, 3, 0 );
Set(buvw, 3, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, 1.0/(sqrt(2) ));
Set(buvy, 0, -1.0/(sqrt(2) ));
Set(buvz, 0, 0 );
Set(buvw, 0, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/(sqrt(2) ));
Set(buvy, 1, -1.0/(sqrt(2) ));
Set(buvz, 1, 0 );
Set(buvw, 1, 0.005716906);


leb_matr.push_back(buffer);
buffer.clear();

//18
//aaa

Set(buvx, 2, 1.0/(sqrt(2) ));
Set(buvy, 2, 1.0/(sqrt(2) ));
Set(buvz, 2, 1.0/(sqrt(2) ));
Set(buvw, 2, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/(sqrt(2) ));
Set(buvy, 3, 1.0/(sqrt(2) ));
Set(buvz, 3, 1.0/(sqrt(2) ));
Set(buvw, 3, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/(sqrt(2) ));
Set(buvy, 0, -1.0/(sqrt(2) ));
Set(buvz, 0, 1.0/(sqrt(2) ));
Set(buvw, 0, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/(sqrt(2) ));
Set(buvy, 1, -1.0/(sqrt(2) ));
Set(buvz, 1, 1.0/(sqrt(2) ));
Set(buvw, 1, 0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 1.0/(sqrt(2) ));
Set(buvy, 2, 1.0/(sqrt(2) ));
Set(buvz, 2, -1.0/(sqrt(2) ));
Set(buvw, 2, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -1.0/(sqrt(2) ));
Set(buvy, 3, 1.0/(sqrt(2) ));
Set(buvz, 3, -1.0/(sqrt(2) ));
Set(buvw, 3, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 1.0/(sqrt(2) ));
Set(buvy, 0, -1.0/(sqrt(2) ));
Set(buvz, 0, -1.0/(sqrt(2) ));
Set(buvw, 0, 0.00557338318);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -1.0/(sqrt(2) ));
Set(buvy, 1, -1.0/(sqrt(2) ));
Set(buvz, 1, -1.0/(sqrt(2) ));
Set(buvw, 1, 0.00557338318);

leb_matr.push_back(buffer);
buffer.clear();

//26
// aab
//second value p1: 0.6713
//q1=p1^2
//lk=sqrt( 1-2 q1*q1)

Set(buvx, 2, 0.6712973);
Set(buvy, 2, 0.6712973);
Set(buvz, 2, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.6712973);
Set(buvy, 3, 0.6712973);
Set(buvz, 3, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.6712973);
Set(buvy, 0, -0.6712973);
Set(buvz, 0, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 0, 0.005608704082);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.6712973);
Set(buvy, 1, -0.6712973);
Set(buvz, 1, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 1, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

//30

Set(buvx, 2, 0.6712973);
Set(buvy, 2, 0.6712973);
Set(buvz, 2, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.6712973);
Set(buvy, 3, 0.6712973);
Set(buvz, 3, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.6712973);
Set(buvy, 0, -0.6712973);
Set(buvz, 0, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 0, 0.005608704082);


leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 1, -0.6712973);
Set(buvy, 1, -0.6712973);
Set(buvz, 1, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvw, 1, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

//34 atras

Set(buvx, 2, 0.6712973);
Set(buvy, 2, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 2, 0.6712973);
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.6712973);
Set(buvy, 3, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 3, 0.6712973);
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.6712973);
Set(buvy, 0, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 0, 0.6712973);
Set(buvw, 0, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 1, -0.6712973);
Set(buvy, 1, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 1, 0.6712973);
Set(buvw, 1, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.6712973);
Set(buvy, 2, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 2, -0.6712973);
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.6712973);
Set(buvy, 3, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 3, -0.6712973);
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 40 atras

Set(buvx, 0, 0.6712973);
Set(buvy, 0, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 0, -0.6712973);
Set(buvw, 0, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.6712973);
Set(buvy, 1, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvz, 1, -0.6712973);
Set(buvw, 1, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 2, 0.6712973);
Set(buvz, 2, 0.6712973);
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 3, 0.6712973);
Set(buvz, 3, 0.6712973);
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 0, -0.6712973);
Set(buvz, 0, 0.6712973);
Set(buvw, 0, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 1, -0.6712973);
Set(buvz, 1, 0.6712973);
Set(buvw, 1, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 2, 0.6712973);
Set(buvz, 2, -0.6712973);
Set(buvw, 2, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 3, 0.6712973);
Set(buvz, 3, -0.6712973);
Set(buvw, 3, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 0, -0.6712973);
Set(buvz, 0, -0.6712973);
Set(buvw, 0, 0.005608704082);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1.0 - 2.0*0.6712973*0.6712973));
Set(buvy, 1, -0.6712973);
Set(buvz, 1, -0.6712973);
Set(buvw, 1, 0.005608704082);


leb_matr.push_back(buffer);
buffer.clear();

//third value  p1^2: 0.289
//q1=p1^2
//lk=sqrt( 1-2 q1*q1) =

// aab(0.00516, 0.289)
//cof a, ob = 0.5158237711805383e-2
// 50 

Set(buvx, 2, 0.289);
Set(buvy, 2, 0.289);
Set(buvz, 2, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.289);
Set(buvy, 3, 0.289);
Set(buvz, 3, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 3, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.289);
Set(buvy, 0, -0.289);
Set(buvz, 0, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 0, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.289);
Set(buvy, 1, -0.289);
Set(buvz, 1, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 1, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.289);
Set(buvy, 2, 0.289);
Set(buvz, 2, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.289);
Set(buvy, 3, 0.289);
Set(buvz, 3, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 3, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.289);
Set(buvy, 0, -0.289);
Set(buvz, 0, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 0, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.289);
Set(buvy, 1, -0.289);
Set(buvz, 1, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvw, 1, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.289);
Set(buvy, 2, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 2, 0.289);
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.289);
Set(buvy, 3, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 3, 0.289);
Set(buvw, 3, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 60 atras

Set(buvx, 0, 0.289);
Set(buvy, 0, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 0, 0.289);
Set(buvw, 0, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.289);
Set(buvy, 1, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 1, 0.289);
Set(buvw, 1, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.289);
Set(buvy, 2, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 2, -0.289);
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.289);
Set(buvy, 3, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 3, -0.289);
Set(buvw, 3, 0.0051582377);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.289);
Set(buvy, 0, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 0, -0.289);
Set(buvw, 0, 0.0051582377);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.289);
Set(buvy, 1, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvz, 1, -0.289);
Set(buvw, 1, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 2, 0.289);
Set(buvz, 2, 0.289);
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 3, 0.289);
Set(buvz, 3, 0.289);
Set(buvw, 3, 0.0051582377);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 0, -0.289);
Set(buvz, 0, 0.289);
Set(buvw, 0, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 1, -0.289);
Set(buvz, 1, 0.289);
Set(buvw, 1, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

// 70 atras

Set(buvx, 2, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 2, 0.289);
Set(buvz, 2, -0.289);
Set(buvw, 2, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 3, 0.289);
Set(buvz, 3, -0.289);
Set(buvw, 3, 0.0051582377);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 0, -0.289);
Set(buvz, 0, -0.289);
Set(buvw, 0, 0.0051582377);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -sqrt(1.0 - 2.0*0.289*0.289));
Set(buvy, 1, -0.289);
Set(buvz, 1, -0.289);
Set(buvw, 1, 0.0051582377);


leb_matr.push_back(buffer);
buffer.clear();

// 74 
//aab ( p1^2 = 0.4447  ,conf = 0.005518)

Set(buvx, 2, 0.4447);
Set(buvy, 2, 0.4447);
Set(buvz, 2, 0.77749);
Set(buvw, 2, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.4447);
Set(buvy, 3, 0.4447);
Set(buvz, 3, 0.77749);
Set(buvw, 3, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.4447);
Set(buvy, 0, -0.4447);
Set(buvz, 0, 0.77749);
Set(buvw, 0, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.4447);
Set(buvy, 1, -0.4447);
Set(buvz, 1, 0.77749);
Set(buvw, 1, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.4447);
Set(buvy, 2, 0.4447);
Set(buvz, 2, -0.77749);
Set(buvw, 2, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.4447);
Set(buvy, 3, 0.4447);
Set(buvz, 3, -0.77749);
Set(buvw, 3, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 80 atras

Set(buvx, 0, 0.4447);
Set(buvy, 0, -0.4447);
Set(buvz, 0, -0.77749);
Set(buvw, 0, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.4447);
Set(buvy, 1, -0.4447);
Set(buvz, 1, -0.77749);
Set(buvw, 1, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.4447);
Set(buvy, 2, 0.77749);
Set(buvz, 2, 0.4447);
Set(buvw, 2, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();


Set(buvx, 3, -0.4447);
Set(buvy, 3, 0.77749);
Set(buvz, 3, 0.4447);
Set(buvw, 3, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.4447);
Set(buvy, 0, -0.77749);
Set(buvz, 0, 0.4447);
Set(buvw, 0, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.4447);
Set(buvy, 1, -0.77749);
Set(buvz, 1, 0.4447);
Set(buvw, 1, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.4447);
Set(buvy, 2, 0.77749);
Set(buvz, 2, -0.4447);
Set(buvw, 2, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.4447);
Set(buvy, 3, 0.77749);
Set(buvz, 3, -0.4447);
Set(buvw, 3, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.4447);
Set(buvy, 0, -0.77749);
Set(buvz, 0, -0.4447);
Set(buvw, 0, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.4447);
Set(buvy, 1, -0.77749);
Set(buvz, 1, -0.4447);
Set(buvw, 1, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

// 90 atras:

Set(buvx, 2, 0.77749);
Set(buvy, 2, 0.4447);
Set(buvz, 2, 0.4447);
Set(buvw, 2, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.77749);
Set(buvy, 3, 0.4447);
Set(buvz, 3, 0.4447);
Set(buvw, 3, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.77749);
Set(buvy, 0, -0.4447);
Set(buvz, 0, 0.4447);
Set(buvw, 0, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.77749);
Set(buvy, 1, -0.4447);
Set(buvz, 1, 0.4447);
Set(buvw, 1, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.77749);
Set(buvy, 2, 0.4447);
Set(buvz, 2, -0.4447);
Set(buvw, 2, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.77749);
Set(buvy, 3, 0.4447);
Set(buvz, 3, -0.4447);
Set(buvw, 3, 0.0055187714);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.77749);
Set(buvy, 0, -0.4447);
Set(buvz, 0, -0.4447);
Set(buvw, 0, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.77749);
Set(buvy, 1, -0.4447);
Set(buvz, 1, -0.4447);
Set(buvw, 1, 0.0055187714);


leb_matr.push_back(buffer);
buffer.clear();

// 98 atras
// aab:
//p1^2 = 0.13, a, C_f, lk o mk = 0.4106777028169394e-2)
//b = sqrt(1-2aa) =   sqrt(1-2 p1^2 p1^2 ) ? = 0.983

Set(buvx, 2, 0.13);
Set(buvy, 2, 0.13);
Set(buvz, 2, 0.983);
Set(buvw, 2, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.13);
Set(buvy, 3, 0.13);
Set(buvz, 3, 0.983);
Set(buvw, 3, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 100 atras

Set(buvx, 0, 0.13);
Set(buvy, 0, -0.13);
Set(buvz, 0, 0.983);
Set(buvw, 0, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.13);
Set(buvy, 1, -0.13);
Set(buvz, 1, 0.983);
Set(buvw, 1, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.13);
Set(buvy, 2, 0.13);
Set(buvz, 2, -0.983);
Set(buvw, 2, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.13);
Set(buvy, 3, 0.13);
Set(buvz, 3, -0.983);
Set(buvw, 3, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.13);
Set(buvy, 0, -0.13);
Set(buvz, 0, -0.983);
Set(buvw, 0, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.13);
Set(buvy, 1, -0.13);
Set(buvz, 1, -0.983);
Set(buvw, 1, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.13);
Set(buvy, 2, 0.983);
Set(buvz, 2, 0.13);
Set(buvw, 2, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.13);
Set(buvy, 3, 0.983);
Set(buvz, 3, 0.13);
Set(buvw, 3, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.13);
Set(buvy, 0, -0.983);
Set(buvz, 0, 0.13);
Set(buvw, 0, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.13);
Set(buvy, 1, -0.983);
Set(buvz, 1, 0.13);
Set(buvw, 1, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

// 110

Set(buvx, 2, 0.13);
Set(buvy, 2, 0.983);
Set(buvz, 2, -0.13);
Set(buvw, 2, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.13);
Set(buvy, 3, 0.983);
Set(buvz, 3, -0.13);
Set(buvw, 3, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.13);
Set(buvy, 0, -0.983);
Set(buvz, 0, -0.13);
Set(buvw, 0, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.13);
Set(buvy, 1, -0.983);
Set(buvz, 1, -0.13);
Set(buvw, 1, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.983);
Set(buvy, 2, 0.13);
Set(buvz, 2, 0.13);
Set(buvw, 2, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.983);
Set(buvy, 3, 0.13);
Set(buvz, 3, 0.13);
Set(buvw, 3, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.983);
Set(buvy, 0, -0.13);
Set(buvz, 0, 0.13);
Set(buvw, 0, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.983);
Set(buvy, 1, -0.13);
Set(buvz, 1, 0.13);
Set(buvw, 1, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.983);
Set(buvy, 2, 0.13);
Set(buvz, 2, -0.13);
Set(buvw, 2, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.983);
Set(buvy, 3, 0.13);
Set(buvz, 3, -0.13);
Set(buvw, 3, 0.004106777);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 120 atras

Set(buvx, 0, 0.983);
Set(buvy, 0, -0.13);
Set(buvz, 0, -0.13);
Set(buvw, 0, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.983);
Set(buvy, 1, -0.13);
Set(buvz, 1, -0.13);
Set(buvw, 1, 0.004106777);

leb_matr.push_back(buffer);
buffer.clear();

// 122 atras
// ab0:
//(0.5051846064614808e-2, 0.3457702197611283e0)
//p1^2 = 0.3457, a, coff = 0.00505)

//b = sqrt(1-2aa) =   sqrt(1- p1^2 p1^2 ) ?
// = sqrt(1 - 0.3457*0.3457) aprox 0.938

Set(buvx, 2, 0.3457);
Set(buvy, 2, 0.938);
Set(buvz, 2, 0);
Set(buvw, 2, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.3457);
Set(buvy, 3, 0.938);
Set(buvz, 3, 0);
Set(buvw, 3, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.3457);
Set(buvy, 0, -0.938);
Set(buvz, 0, 0);
Set(buvw, 0, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.3457);
Set(buvy, 1, -0.938);
Set(buvz, 1, 0);
Set(buvw, 1, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.938);
Set(buvy, 2, 0.3457);
Set(buvz, 2, 0);
Set(buvw, 2, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.938);
Set(buvy, 3, 0.3457);
Set(buvz, 3, 0);
Set(buvw, 3, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.938);
Set(buvy, 0, -0.3457);
Set(buvz, 0, 0);
Set(buvw, 0, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.938);
Set(buvy, 1, -0.3457);
Set(buvz, 1, 0);
Set(buvw, 1, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

// 130 atras

Set(buvx, 2, 0.34577);
Set(buvy, 2, 0);
Set(buvz, 2, 0.938);
Set(buvw, 2, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.34577);
Set(buvy, 3, 0);
Set(buvz, 3, 0.938);
Set(buvw, 3, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.34577);
Set(buvy, 0, 0);
Set(buvz, 0, -0.938);
Set(buvw, 0, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.34577);
Set(buvy, 1, 0);
Set(buvz, 1, -0.938);
Set(buvw, 1, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.938);
Set(buvy, 2, 0);
Set(buvz, 2, 0.34577);
Set(buvw, 2, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.938);
Set(buvy, 3, 0);
Set(buvz, 3, 0.34577);
Set(buvw, 3, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.938);
Set(buvy, 0, 0);
Set(buvz, 0, -0.34577);
Set(buvw, 0, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.938);
Set(buvy, 1, 0);
Set(buvz, 1, -0.34577);
Set(buvw, 1, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0.34577);
Set(buvz, 2, 0.938);
Set(buvw, 2, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -0.34577);
Set(buvz, 3, 0.938);
Set(buvw, 3, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 140 atras

Set(buvx, 0, 0);
Set(buvy, 0, 0.34577);
Set(buvz, 0, -0.938);
Set(buvw, 0, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, -0.34577);
Set(buvz, 1, -0.938);
Set(buvw, 1, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0.938);
Set(buvz, 2, 0.34577);
Set(buvw, 2, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, 0);
Set(buvy, 3, -0.938);
Set(buvz, 3, 0.34577);
Set(buvw, 3, 0.005051846);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0);
Set(buvy, 0, 0.938);
Set(buvz, 0, -0.34577);
Set(buvw, 0, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, 0);
Set(buvy, 1, -0.938);
Set(buvz, 1, -0.34577);
Set(buvw, 1, 0.005051846);

leb_matr.push_back(buffer);
buffer.clear();

//abc:
//0.5530248916233094e-2, 0.1590417105383530e0, 0.8360360154824589e0)
//p1^2 (1), lk(1) o mk(1) segun = 0.16, a, coff = 0.00553)
//p2^2 (2), lk(2) o mk(2) segun  = 0.836
//c = sqrt(1.0 - a*a - b*b) = sqrt(1.0 - 0.16*0.16 - 0.836*0.836 ) aprox 0.525
// 146 atras

Set(buvx, 2, 0.16);
Set(buvy, 2, 0.836);
Set(buvz, 2, 0.52512);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.16);
Set(buvy, 3, 0.836);
Set(buvz, 3, 0.52512);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.16);
Set(buvy, 0, -0.836);
Set(buvz, 0, 0.52512);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.16);
Set(buvy, 1, -0.836);
Set(buvz, 1, 0.52512);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

// 150 atras

Set(buvx, 2, 0.16);
Set(buvy, 2, 0.836);
Set(buvz, 2, -0.52512);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.16);
Set(buvy, 3, 0.836);
Set(buvz, 3, -0.52512);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.16);
Set(buvy, 0, -0.836);
Set(buvz, 0, -0.52512);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.16);
Set(buvy, 1, -0.836);
Set(buvz, 1, -0.52512);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.16);
Set(buvy, 2, 0.52512);
Set(buvz, 2, 0.836);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.16);
Set(buvy, 3, 0.52512);
Set(buvz, 3, 0.836);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.16);
Set(buvy, 0, -0.52512);
Set(buvz, 0, 0.836);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.16);
Set(buvy, 1, -0.52512);
Set(buvz, 1, 0.836);
Set(buvw, 1, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.16);
Set(buvy, 2, 0.52512);
Set(buvz, 2, -0.836);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.16);
Set(buvy, 3, 0.52512);
Set(buvz, 3, -0.836);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 160 atras

Set(buvx, 0, 0.16);
Set(buvy, 0, -0.52512);
Set(buvz, 0, -0.836);
Set(buvw, 0, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.16);
Set(buvy, 1, -0.52512);
Set(buvz, 1, -0.836);
Set(buvw, 1, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.836);
Set(buvy, 2, 0.16);
Set(buvz, 2, 0.52512);
Set(buvw, 2, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.836);
Set(buvy, 3, 0.16);
Set(buvz, 3, 0.52512);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.836);
Set(buvy, 0, -0.16);
Set(buvz, 0, 0.52512);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.836);
Set(buvy, 1, -0.16);
Set(buvz, 1, 0.52512);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.836);
Set(buvy, 2, 0.16);
Set(buvz, 2, -0.52512);
Set(buvw, 2, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.836);
Set(buvy, 3, 0.16);
Set(buvz, 3, -0.52512);
Set(buvw, 3, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.836);
Set(buvy, 0, -0.16);
Set(buvz, 0, -0.52512);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.836);
Set(buvy, 1, -0.16);
Set(buvz, 1, -0.52512);
Set(buvw, 1, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

// 170 atras 

Set(buvx, 2, 0.836);
Set(buvy, 2, 0.52512);
Set(buvz, 2, 0.16);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.836);
Set(buvy, 3, 0.52512);
Set(buvz, 3, 0.16);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.836);
Set(buvy, 0, -0.52512);
Set(buvz, 0, 0.16);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.836);
Set(buvy, 1, -0.52512);
Set(buvz, 1, 0.16);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.836);
Set(buvy, 2, 0.52512);
Set(buvz, 2, -0.16);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.836);
Set(buvy, 3, 0.52512);
Set(buvz, 3, -0.16);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.836);
Set(buvy, 0, -0.52512);
Set(buvz, 0, -0.16);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.836);
Set(buvy, 1, -0.52512);
Set(buvz, 1, -0.16);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.52512);
Set(buvy, 2, 0.16);
Set(buvz, 2, 0.836);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.52512);
Set(buvy, 3, 0.16);
Set(buvz, 3, 0.836);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

// 180

Set(buvx, 0, 0.52512);
Set(buvy, 0, -0.16);
Set(buvz, 0, 0.836);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.52512);
Set(buvy, 1, -0.16);
Set(buvz, 1, 0.836);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.52512);
Set(buvy, 2, 0.16);
Set(buvz, 2, -0.836);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.52512);
Set(buvy, 3, 0.16);
Set(buvz, 3, -0.836);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.52512);
Set(buvy, 0, -0.16);
Set(buvz, 0, -0.836);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.52512);
Set(buvy, 1, -0.16);
Set(buvz, 1, -0.836);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0.52512);
Set(buvy, 2, 0.836);
Set(buvz, 2, 0.16);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.52512);
Set(buvy, 3, 0.836);
Set(buvz, 3, 0.16);
Set(buvw, 3, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

Set(buvx, 0, 0.52512);
Set(buvy, 0, -0.836);
Set(buvz, 0, 0.16);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.52512);
Set(buvy, 1, -0.836);
Set(buvz, 1, 0.16);
Set(buvw, 1, 0.005530249);


leb_matr.push_back(buffer);
buffer.clear();

// 190 atras

Set(buvx, 2, 0.52512);
Set(buvy, 2, 0.836);
Set(buvz, 2, -0.16);
Set(buvw, 2, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 3, -0.52512);
Set(buvy, 3, 0.836);
Set(buvz, 3, -0.16);
Set(buvw, 3, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();


Set(buvx, 0, 0.52512);
Set(buvy, 0, -0.836);
Set(buvz, 0, -0.16);
Set(buvw, 0, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 1, -0.52512);
Set(buvy, 1, -0.836);
Set(buvz, 1, -0.16);
Set(buvw, 1, 0.005530249);

leb_matr.push_back(buffer);
buffer.clear();

Set(buvx, 2, 0);
Set(buvy, 2, 0);
Set(buvz, 2, 0);
Set(buvw, 2, 0);

Set(buvx, 3, 0);
Set(buvy, 3, 0);
Set(buvz, 3, 0);
Set(buvw, 3, 0);

lebedev_buffer.push_back(buvx);
lebedev_buffer.push_back(buvy);
lebedev_buffer.push_back(buvz);
lebedev_buffer.push_back(buvw);

lebedev_total.push_back(lebedev_buffer);
lebedev_buffer.clear();

break;
}

return lebedev_total;
}

VECCORE_ATT_HOST_DEVICE
Double_v Power(Double_v base, G4int ex){
//until 10 int exponent
if( ex < 0 ){

return 1.0/Power( base, abs( ex) );

}else{   


MaskI_v p0 (ex == 0 );
if(MaskFull(p0)) {
return 1.0;
}else{
MaskI_v p1 (ex == 1 );

if(MaskFull(p1)) {
return base;
}else{

MaskI_v p2 (ex == 2 );
if(MaskFull(p2)) {
return base*base;
}else{

MaskI_v p3 (ex == 3 );
if(MaskFull(p3)) {
return base*base*base;
}else{
MaskI_v p4 (ex == 4 );
if(MaskFull(p4)) {
return base*base*base*base;
}else{

MaskI_v p5 (ex == 5 );
if(MaskFull(p5)) {
return base*base*base*base*base;
}else{

MaskI_v p6 (ex == 6 );
if(MaskFull(p6)) {
return base*base*base*base*base*base;
}else{

MaskI_v p7 (ex == 7 );
if(MaskFull(p7)) {
return base*base*base*base*base*base*base;
}else{

MaskI_v p8 (ex == 8 );
if(MaskFull(p8)) {
return base*base*base*base*base*base*base*base;


}else{
MaskI_v p9 (ex == 9 );
if(MaskFull(p9 )) {
return base*base*base*base*base*base*base*base*base;

}else{
MaskI_v p10 (ex == 10 );
if(MaskFull(p10 )) {
return base*base*base*base*base*base*base*base*base*base;

}else{

MaskI_v paro (ex > 10 );
if(MaskFull(p10 )) {
printf("++ naive power \n");
return 1.0;
}


}

}

}
  
}

}

}

}

}

}  

}

}

//check negative:
}

}


Double_v Correlation_Interchange::evaluation_functionV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);

return N*( Power( X - X_At, L[0] ))*(Power( Y - Y_At, L[1] ))*
(Power( Z - Z_At, L[2]) )*vecCore::math::Exp(-Alpha*rV);
}


Double_v Correlation_Interchange::evaluation_partial_functionxV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );


if(   MaskFull( dif_x )   ) {

X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;
}

}

}


Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( L[0]*(Power( X - X_At , L[0] - 1 ))*
(Power( Y - Y_At, L[1] ))*(Power( Z - Z_At , L[2] ) ) - 2.0*Alpha*(Power(X - X_At, L[0] + 1  ))*
(Power( Y - Y_At , L[1] ))*(Power( Z - Z_At , L[2] ) )     	)  );

MaskD_v nan_v ( Ev  == true  );
if(   MaskFull( nan_v )   ) {
printf("vxxx \n");
}

return Ev;

}

Double_v Correlation_Interchange::evaluation_partial_functionyV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_y ( (Y - Y_At)  == 0  );

if(   MaskFull( dif_y )   ) {

Y = Y + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( L[1]*(Power( X - X_At, L[0] ))*(Power( Y - Y_At , L[1] - 1  ))*
(Power( Z - Z_At , L[2] )) - 2.0*Alpha*(Power( X - X_At, L[0] ))*(Power( Y - Y_At , L[1] + 1 ))*(Power( Z - Z_At , L[2] ))  )  ) ;

MaskD_v nan_v ( Ev  == true  );

return Ev;
}



Double_v Correlation_Interchange::evaluation_partial_functionzV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_z )   ) {
Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

//Get(C , k ), Get(VB[id] , k ) , k,  contv);
if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;
//x = x + 0.000001; 
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( L[2]*(Power( X - X_At ,L[0] ))*(Power( Y - Y_At, L[1] ))*
(Power( Z - Z_At, L[2] - 1 )) - 2.0*Alpha*(Power( X - X_At, L[0] ))*(Power( Y - Y_At, L[1] ))*(Power ( Z - Z_At, L[2] + 1 ) ) )  );

MaskD_v nan_v ( Ev  == true  );
if(   MaskFull( nan_v )   ) {
printf("non-nand \n");
}

return Ev;
}


Double_v Correlation_Interchange::evaluation_double_partial_function_xV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );

if(   MaskFull( dif_x )   ) {
X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
//printf("WA %f C %f VB[id] %f en k %d contv %d \n",  Get(WA , k ) ,
//Get(C , k ), Get(VB[id] , k ) , k,  contv);

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;
//x = x + 0.000001; 
}

}

}


Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*(Power( Y - Y_At, L[1]) )*(Power( Z - Z_At, L[2] ))*
(  ( L[0]*L[0] - L[0] )*( Power( X - X_At, L[0] - 2 ) ) - 2.0*Alpha*(Power( X - X_At, L[0] ))*
(2*L[0] + 1) + 4.0*Alpha*Alpha*(Power( X - X_At, L[0] + 2 ))    )  );

MaskD_v nan_v ( Ev  == true  );
//if(   MaskFull( nan_v )   ) {
//printf("xxxx \n");
//}

return Ev;
}


Double_v Correlation_Interchange::evaluation_double_partial_function_yV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_y ( (Y - Y_At)  == 0  );

if(   MaskFull( dif_y )   ) {
Y = Y + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*(Power( X - X_At, L[0]) )*(Power( Z - Z_At, L[2] ))*
(  ( L[1]*L[1] - L[1] )*(Power( Y - Y_At, L[1] - 2  )) - 2.0*Alpha*(Power( Y - Y_At, L[1] ))*
(2*L[1] + 1) + 4.0*Alpha*Alpha*(Power ( Y - Y_At, L[1] + 2 ))  )   );

MaskD_v nan_v ( Ev  == true  );

return Ev;
}

Double_v Correlation_Interchange::evaluation_double_partial_function_zV ( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_z )   ) {
Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;

}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*(Power( X - X_At, L[0]))*(Power( Y - Y_At, L[1] ))*
( ( L[2]*L[2] - L[2] )*(Power ( Z - Z_At , L[2] - 2 )) -2.0*Alpha*(Power( Z - Z_At , L[2]  ))*
(2*L[2] + 1) + 4*Alpha*Alpha*(Power( Z - Z_At,  L[2] + 2  ))  )  );

MaskD_v nan_v ( Ev  == true  );

return Ev;
}



Double_v Correlation_Interchange::evaluation_partial_function_dydxV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){   

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );
MaskD_v dif_y ( (Y - Y_At)  == 0  );

if(   MaskFull( dif_x )   ) {

X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;

}

}

}


if(   MaskFull( dif_y )   ) {

Y = Y + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;

}

}

}


Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( Power( Z - Z_At , L[2] ))*
( L[0]*(Power( X - X_At, L[0] - 1)) -2.0*Alpha*(Power( X - X_At, L[0] + 1 ))  )*
( L[1]*(Power( Y - Y_At, L[1] - 1 )) -2.0*Alpha*(Power( Y - Y_At, L[1] + 1 ) )  )  );

MaskD_v nan_v ( Ev  == true  );

return Ev;
}


Double_v Correlation_Interchange::evaluation_partial_function_dzdxV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){   

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_x )   ) {
X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;
}

}

}


if(   MaskFull( dif_z )   ) {
Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( Power( Y - Y_At, L[1] ) )*
( L[0]*(Power( X - X_At, L[0] - 1 )) - 2.0*Alpha*(Power ( X - X_At, L[0] + 1 ))  )*
( L[2]*(Power( Z - Z_At, L[2] - 1 )) - 2.0*Alpha*(Power ( Z - Z_At, L[2] + 1 ))  ) );   

MaskD_v nan_v ( Ev  == true  );

return Ev;
}


Double_v Correlation_Interchange::evaluation_partial_function_dxdyV( G4double X_At,
G4double Y_At, G4double Z_At, G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){  

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );
MaskD_v dif_y ( (Y - Y_At)  == 0  );

if(   MaskFull( dif_x )   ) {
X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;
}

}

}

if(   MaskFull( dif_y )   ) {

Y = Y + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*(  L[1]*(Power( Y - Y_At, L[1] -1 ) ) -
2.0*Alpha*(Power( Y - Y_At, L[1] + 1 )) )*(Power( Z - Z_At, L[2] ))*
( L[0]*(Power( X - X_At , L[0] -1 )) - 2.0*Alpha*(Power( X - X_At , L[0] + 1 ))   ) );

MaskD_v nan_v ( Ev  == true  );

return Ev;
}

Double_v Correlation_Interchange::evaluation_partial_function_dxdzV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){  


Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_x ( (X - X_At)  == 0  );
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_x )   ) {

X = X + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (X[k] - X_At) == 0  ){
Set( X , k , X[k] + 0.000001 ) ;
}

}

}


if(   MaskFull( dif_z )   ) {

Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( Power( Y - Y_At, L[1] ))*
( L[2]*(Power( Z - Z_At, L[2] - 1 )) - 2.0*Alpha*(Power( Z - Z_At,  L[2] + 1  )) )*
( L[0]*(Power( X - X_At, L[0] - 1 )) - 2.0*Alpha*(Power( X - X_At, L[0] + 1 )) ) ) ;

MaskD_v nan_v ( Ev  == true  );

return Ev;
}


Double_v Correlation_Interchange::evaluation_partial_function_dzdyV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){  

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_y ( (Y - Y_At)  == 0  );
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_y )   ) {
Y = Y + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;

}

}

}

if(   MaskFull( dif_z )   ) {

Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( Power( X - X_At, L[0] ) )*
( L[1]*(Power( Y - Y_At, L[1] - 1  )) - 2.0*Alpha*(Power( Y - Y_At, L[1] + 1 ))  )*
( L[2]*(Power( Z - Z_At, L[2] - 1  )) - 2.0*Alpha*(Power( Z - Z_At, L[2] + 1 )) ) );             

MaskD_v nan_v ( Ev  == true  );

return Ev;
}


Double_v Correlation_Interchange::evaluation_partial_function_dydzV( G4double X_At,
G4double Y_At, G4double Z_At,  G4double Alpha, G4double N,
vector<G4int> L, Double_v X, Double_v Y, Double_v Z ){  

Double_v rV = (X - X_At)*(X - X_At) + (Y - Y_At)*( Y - Y_At) + (Z - Z_At)*(Z - Z_At);
MaskD_v dif_y ( (Y - Y_At)  == 0  );
MaskD_v dif_z ( (Z - Z_At)  == 0  );

if(   MaskFull( dif_y )   ) {
Y = Y + 0.000001; 
}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {

if( (Y[k] - Y_At) == 0  ){
Set( Y , k , Y[k] + 0.000001 ) ;

}

}

}

if(   MaskFull( dif_z )   ) {

Z = Z + 0.000001; 

}else{

for (G4int k = 0; k < kVecLenD ; ++k ) {
if( (Z[k] - Z_At) == 0  ){
Set( Z , k , Z[k] + 0.000001 ) ;
}

}

}

Double_v Ev = N*( (vecCore::math::Exp(-Alpha*rV))*( Power( X - X_At, L[0] ) )*
( L[2]*(Power( Z - Z_At, L[2] - 1 )) - 2.0*Alpha*(Power( Z - Z_At, L[2] + 1 )) )*
( L[1]*(Power( Y - Y_At, L[1] - 1 )) - 2.0*Alpha*(Power( Y - Y_At,  L[1] + 1 )) )  );

MaskD_v nan_v ( Ev  == true  );

return Ev;
}



vector<Double_v> Correlation_Interchange::becke_eV( Atom& atomo_puntos ,  Molecule& Mol, G4int& cont,
G4int id , G4int functional ){

vector<Residue *>::iterator res_i;
vector<Atom *>::iterator atr_i;

vector<Residue *>::iterator res_j;
vector<Atom *>::iterator atr_j;

vector<G4double> pp;
vector<Double_v> VB;
vector<Double_v> BeckV;
vector<G4double> pesos;

vector<Residue *>::iterator res_k;
vector<Atom *>::iterator atr_k;

Double_v X,Y,Z;
G4double R, uab, a;
Double_v  B, Bx, By, Bz, RA, RB, RV, M, F, VAR, C, WA;

Double_v Ev = Double_v (0.0);
Double_v Evx = Double_v (0.0);
Double_v Evy = Double_v (0.0);
Double_v Evz = Double_v (0.0);

//G4int cont = 0;
G4int dd = 0;
G4int vv;
G4int size_puntos = atomo_puntos.atomic_grid.size();

Double_v DxDx, DyDy, DzDz;
Double_v DxDy, DxDz;
Double_v DyDx, DyDz;
Double_v DzDx, DzDy;

Double_v SDxDx = Double_v (0.0);
Double_v SDyDy = Double_v (0.0);
Double_v SDzDz = Double_v (0.0);
Double_v SDxDy = Double_v (0.0);
Double_v SDxDz = Double_v (0.0);
Double_v SDyDx = Double_v (0.0);
Double_v SDyDz = Double_v (0.0);
Double_v SDzDx = Double_v (0.0);
Double_v SDzDy = Double_v (0.0);

//for through EMl vector
for(int i = 0; i < size_puntos ; i++) {

G4int size_mallaV = atomo_puntos.atomic_grid[i].lebV.size();

for(int j = 0; j < size_mallaV ; j++) {

X = atomo_puntos.atomic_grid[i].lebV[j][0];
Y = atomo_puntos.atomic_grid[i].lebV[j][1];
Z = atomo_puntos.atomic_grid[i].lebV[j][2];

Double_v XX_A;
Double_v YY_A;
Double_v ZZ_A;

Double_v XX_B;
Double_v YY_B;
Double_v ZZ_B;

G4int atm_cont = 0;


for(res_i = Mol.Residuos_cadena.begin();
    res_i != Mol.Residuos_cadena.end() ; ++res_i ){

for(atr_i = (*res_i)->Lista_atoms.begin();
     atr_i != (*res_i)->Lista_atoms.end(); ++ atr_i ){

XX_A = Power( (*atr_i)->fX - X  ,2);
YY_A = Power( (*atr_i)->fY - Y  ,2);
ZZ_A = Power( (*atr_i)->fZ - Z  ,2);

RA =  vecCore::math::Sqrt( XX_A + YY_A + ZZ_A );

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

XX_B = Power( (*atr_j)->fX - X  ,2);
YY_B = Power( (*atr_j)->fY - Y  ,2);
ZZ_B = Power( (*atr_j)->fZ - Z  ,2);


RB =  vecCore::math::Sqrt( XX_B + YY_B + ZZ_B );

//////////////////////////////////////////////////////////////////////////////

R = sqrt(  ( (*atr_i)->fX - (*atr_j)->fX)*( (*atr_i)->fX - (*atr_j)->fX) +
( (*atr_i)->fY - (*atr_j)->fY)*( (*atr_i)->fY - (*atr_j)->fY) +
( (*atr_i)->fZ - (*atr_j)->fZ)*( (*atr_i)->fZ - (*atr_j)->fZ) );

M = (RA - RB)/(R);

// if heteronuclear
if(  (*atr_i)->fZat  !=   (*atr_j)->fZat  ){

uab = ( (*atr_i)->Bragg_radius -
(*atr_j)->Bragg_radius )/((*atr_i)->Bragg_radius +
(*atr_j)->Bragg_radius );

a = uab/((uab*uab)-1);
//  0.5 prevent div /0 AS 0.45 
F = 0.0 ;

if(a  > 0.5  ){
a =  0.5;
}

if(a  < -0.5  ){
a = -0.5;
}

F = M + a*(1-M*M);
// if homonuclear
}else{
F = M;
}

for(int ii = 0; ii < 3 ; ii++) {
F = F*(3.0 - F*F)/2.0;
}

VAR = VAR*(1 - F )/2.0;
}//fin if

}
} // closing atoms in j

VB.push_back(VAR);
VAR = 1.0;

}
} // closing  atom i

//closing ij atom
C = 0.0;
vv = VB.size();
for(int p = 0; p < vv ; p++) {
C = C + VB[p];
}

//prevent anan
////////////////////////////////////////////////////////
MaskD_v Acc ( C == 0 );

if(   MaskFull( Acc )   ) {
WA  = 0.0;
}else{


if(  !MaskFull( Acc ) ){
WA = VB[id] / C;

}else{


}

for (G4int k = 0; k < kVecLenD ; ++k ) {

if(C[k]  == 0   ){
Set( WA, k , 0);
}

}

}


BeckV.push_back( WA*atomo_puntos.atomic_grid[i].lebV[j][3] ); 

VB.clear();
G4int cont_int = 0;

for(res_k = Mol.Residuos_cadena.begin();
    res_k != Mol.Residuos_cadena.end() ; ++res_k){

 for(atr_k = (*res_k)->Lista_atoms.begin();
     atr_k != (*res_k)->Lista_atoms.end(); ++ atr_k){

  for((*atr_k)->shell1 = (*atr_k)->orbitals.begin();
      (*atr_k)->shell1 != (*atr_k)->orbitals.end(); ++ (*atr_k)->shell1){

for (int k = 0; k < 3; k++) {

// (*(*atr_k)->shell1).NAC[k].A, (*(*atr_k)->shell1).NAC[k].N,
B = evaluation_functionV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

Ev = Ev + B*( Get(  (*(*atr_k)->shell1).orbitalsConst , k ) );

if(functional != 0 ){

Bx = evaluation_partial_functionxV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );


MaskD_v nan_bx ( Bx  == true  );

Evx = Evx + Bx*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

By = evaluation_partial_functionyV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

Evy = Evy + By*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

Bz = evaluation_partial_functionzV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

Evz = Evz + Bz*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

///////////////////////////////////////////////////////////////////////////
// double derivatives

DxDx = evaluation_double_partial_function_xV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

MaskD_v nan_dxdx ( DxDx  == true  );
if(   MaskFull( nan_dxdx )   ) {
printf("valio \n");
}

SDxDx = SDxDx + DxDx*( Get( (*(*atr_k)->shell1).orbitalsConst , k )  );


DyDy = evaluation_double_partial_function_yV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDyDy = SDyDy + DyDy*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DzDz = evaluation_double_partial_function_zV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDzDz = SDzDz + DzDz*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

//

DxDy = evaluation_partial_function_dxdyV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDxDy = SDxDy + DxDy*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DxDz = evaluation_partial_function_dxdzV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDxDz = SDxDz + DxDz*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DyDx = evaluation_partial_function_dydxV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDyDx = SDyDx + DyDx*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DyDz = evaluation_partial_function_dydzV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

MaskD_v nan_dydz ( DyDz  == true  );

SDyDz = SDyDz + DyDz*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DzDx = evaluation_partial_function_dzdxV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDzDx = SDzDx + DzDx*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

DzDy = evaluation_partial_function_dzdyV( (*atr_k)->fX, (*atr_k)->fY,  (*atr_k)->fZ,
Get(  (*(*atr_k)->shell1).orbitalsAlpha , k) , Get( (*(*atr_k)->shell1).orbitalsNormals , k ),
(*(*atr_k)->shell1).l, X ,Y, Z );

SDzDy = SDzDy + DzDy*( Get( (*(*atr_k)->shell1).orbitalsConst , k ) );

}


}


///////////////////////////////////////////////////////////////////////////

//for (G4int k = 0; k < kVecLenD - 1; ++k ) {
//ev = ev + Get(B, k)*Get(    (*(*atr_k)->shell1).orbitalsConst  , k );
//}

////////////////////////////////////////////////////////////////////////

(*(*atr_k)->shell1).EvalV.push_back(Ev);

if(functional != 0 ){

(*(*atr_k)->shell1).X_G.push_back(Evx);
(*(*atr_k)->shell1).Y_G.push_back(Evy);
(*(*atr_k)->shell1).Z_G.push_back(Evz);

(*(*atr_k)->shell1).DXX.push_back(SDxDx);
(*(*atr_k)->shell1).DYY.push_back(SDyDy);
(*(*atr_k)->shell1).DZZ.push_back(SDzDz);

(*(*atr_k)->shell1).DXY.push_back(SDxDy);
(*(*atr_k)->shell1).DXZ.push_back(SDxDz);

(*(*atr_k)->shell1).DYX.push_back(SDyDx);
(*(*atr_k)->shell1).DYZ.push_back(SDyDz);

(*(*atr_k)->shell1).DZX.push_back(SDzDx);
(*(*atr_k)->shell1).DZY.push_back(SDzDy);

}


Ev = Double_v (0.0);
Evx = Double_v (0.0);
Evy = Double_v (0.0);
Evz = Double_v (0.0);

SDxDx = Double_v (0.0);
SDyDy = Double_v (0.0);
SDzDz = Double_v (0.0);

SDxDy = Double_v (0.0);
SDxDz = Double_v (0.0);

SDyDx = Double_v (0.0);
SDyDz = Double_v (0.0);

SDzDx = Double_v (0.0);
SDzDy = Double_v (0.0);


cont_int++;
}
}
} 

cont++;
}// j
}// i


return BeckV;
}


void Correlation_Interchange::Euler_Mac_Lebedev_SG1 ( Atom& atomo){

vector<vector<Double_v>>::iterator ColV;
vector<Double_v>::iterator RowV;

G4int Nr = 50;
G4int siz_puntos;
G4double wir =0;
G4double ri = 0;
G4double bf = 0;
G4int No = 0;
G4double RE = 0;
//G4int cont = 0;
//G4int contv= 0;

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
// 

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

if( ( atomo.Bohr_radius == 1.7862) || ( atomo.Bohr_radius==  1.6489)  ){

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

// default
No = 86;
}
}
}
}

}
}
}

points.Nomega = No;
points.lebV  = lebedev_weightsV (No);
siz_puntos = points.lebV.size();

for(int i = 0; i <  siz_puntos ; i++) {

//x is  0:   [i][0]
points.lebV[i][0] = points.lebV[i][0]*ri + atomo.fX;
points.lebV[i][1] = points.lebV[i][1]*ri + atomo.fY;
points.lebV[i][2] = points.lebV[i][2]*ri + atomo.fZ;
points.lebV[i][3] = points.lebV[i][3]*wir*(4.0)*pi;

}

atomo.atomic_grid.push_back(points);

}

return;
}




vector<Double_v>  Correlation_Interchange::GridV( Molecule& Mol ,  G4int funcional ){

vector<Atom *>::iterator atr_main;
vector<Atom *>::iterator atr_row;
vector<Residue *>::iterator res_main;
vector<Residue *>::iterator res_row;
vector<Double_v>  Wesos;
vector<Double_v>  WVB;

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

WVB =  becke_eV (  (*(*atr_main)) ,  Mol , sis, cont, funcional );

Mol.No_points = Mol.No_points + sis;

Wesos.insert(Wesos.end(), WVB.begin(), WVB.end());

cont++;
}
}

return Wesos;
}
