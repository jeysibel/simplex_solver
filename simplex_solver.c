/*
 * Projeto da Disciplina Pesquisa Operacional
 *  professor: Anand Subramanian
 *  @author Jeysibel de Sousa Dantas
 *
 * Simplex Solver
 *    Solver do Método Simplex
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <float.h>

#define VARNAME_SIZE 10

#ifdef DEBUG
# define DEBUG_PRINT(x) printf x
#else
# define DEBUG_PRINT(x) do {} while (0)
#endif

//armazena o numero de variáveis inicial do problema
int qtyVariables;
//armazena o numero de variaveis de folga necessarias a normalização do problema
int qtySlackVariables;
//armazena o numero de váriaveis necessarias a normalização do problema
int qtyArtificialVariables;

//armazena a quantidade de regras de restricao do problema, excluindo as de positividade
int qtyRestrictions;
//valor total acrescido das regras de não negatividade
int qtyRestrictionsTotal;

//armazena o numero total de colunas do tableau ( Xn + XnF + XnA + b) 
int qtyCollumnsTableau;
//armazena o numero total de linhas do tableau (FO + (Restricoes exceto de +) )
int qtyLinesTableau;

//armazena o numero total de linhas to tableau incluindo as resticoes de não negatividade
int qtyLinesRestrictionsTable;

//armazena os nomes das variaveis para cada coluna do tableau
char **collumnNames;
//armazena os nomes das variaveis de saida para cada linha do tableau
char **lineNames;

//Função objetivo normalizada
double *rawObjectiveFunction;
//Restrições Básicas normalizadas
double **rawRestrictions;

//Função objetivo normalizada
double *normObjectiveFunction;
//Restrições Básicas normalizadas
double **normRestrictions;
//Tableau composto de todas as restrições e FO normalizadas para o simplex
double **tableau;


int iterations = 0;

int max_iterations = 0;



/*
*  esta função é responsável por imprimir em tela os vários status do tableau durante iterações
*/
void print_tableau_state(){
  
  printf("\n ===Tableau Status after %i iterations===\t\t\t\n\n| Ln\\Xn  |\t", iterations);

  for (int i=0; i< (qtyCollumnsTableau); i++){
    char *collumnItem = *(collumnNames + i);
    //DEBUG_PRINT(( "\nCollumNames addr: %p , CollumnItem: %p , i: %i \n", collumnNames, collumnItem, i));
    printf("%s\t|\t", collumnItem);
  }
  printf("\n");

  for (int i=0; i < qtyLinesTableau ;i++){
    char *lineName = *(lineNames + i);
    printf("| %s\t |", lineName );
    double *tableauLine = *(tableau + i);
    //DEBUG_PRINT(( "\nLineName addr: %p , tableauLine: %p , i: %i \n", lineName, tableauLine, i));
    for (int j=0; j < qtyCollumnsTableau; j++ ){
        printf(" %f\t|", tableauLine[j] );
    }
    printf("\n");
  }
}


void print_tableau_variables(){
  
  printf("\n ===Tableau Variables after %i iterations: ===\t\t\t\n\n\t", iterations);

  //qtyCollumns represents the max possible size
  double *basicVarResults = malloc (sizeof (double)  * qtyCollumnsTableau );
  char **basicVarNames = malloc (sizeof (char*) * qtyCollumnsTableau);
  char **nonBasicVarNames = malloc (sizeof (char*) * qtyCollumnsTableau);

  int basicVarIndex = 0;
  int nonBasicVarIndex = 0;

  int searchRange = qtyCollumnsTableau -1;
  for (int i=0; i< searchRange; i++){
    bool isBasicVariable = false;
    char *collumnItem = *(collumnNames + i);
    DEBUG_PRINT(("\n\n >>>>Checking %s \n", collumnItem));
    double bLineResult = -1.0;
    for (int j=0; j < qtyLinesTableau ;j++){
      double *lineValue = *(tableau + j);        
      double currentValue = *(lineValue + i);
       DEBUG_PRINT(("for current Value: %f at line %i\n", currentValue , j));
      if (currentValue != 1.0 && currentValue != 0.0){
        DEBUG_PRINT(("\n%s Discarded as Basic Variable: %f\n", collumnItem ,  currentValue));
        isBasicVariable = false;
        break;
      }else if( currentValue == 1.0){
        DEBUG_PRINT(("\n%s Temporarily Accepted as Basic Variable: %f\n", collumnItem, currentValue));
        isBasicVariable = true;

        DEBUG_PRINT(("Adding %f as result for %s\n", *(lineValue +  searchRange), collumnItem ));
        bLineResult = *(lineValue +  searchRange);
      }
    
    }

    if( isBasicVariable){
      *(basicVarNames + basicVarIndex) = collumnItem;
      *(basicVarResults + basicVarIndex) = bLineResult;
      basicVarIndex++;  
    } else{
      *(nonBasicVarNames + nonBasicVarIndex) = collumnItem;
      nonBasicVarIndex++;
    }

  }

  printf("\n\t * Estimated Z value: *\n\n\t\t %s = %f\n",  *(basicVarNames), *(basicVarResults) );
  printf("\n\t * Basic Variables *\n");
  for (int i = 1; i < basicVarIndex; i++ ){
    printf("\n\t\t %s = %f\n", *(basicVarNames+i), *(basicVarResults+i) );
  }
  printf("\n\t * Non Basic Variables *\n");
  for (int i = 0; i < nonBasicVarIndex ; i++ ){
    printf("\n\t\t %s = 0.0\n", *(nonBasicVarNames+i) );
  }


}


void print_tableau_equations(){

  int qtyTerms =  qtyVariables + 1;
  int qtyVar = qtyVariables;
  int qtyRests =  qtyRestrictions;

  printf("\n ===Init Problem Solver For :=== \n\n");

  //DEBUG_PRINT(("\n qtyTerms: %d  qtyVar: %d qtyRests: %d \n", qtyTerms, qtyVar, qtyRests));

  for ( int i = 0; i < qtyTerms; i++) {
    //DEBUG_PRINT(("\n i: %d  rawObjectiveFunction + i: %f String: ", i, *(rawObjectiveFunction + i) ));

    if (i == 0 ) {
      printf("   Z = ");
    } else if ( i < qtyTerms - 1 ){
      printf(" %.0fX%i", *(rawObjectiveFunction + i), i );
      printf(" + ");
    } else {
      printf(" %.0fX%i ", *(rawObjectiveFunction + i), i );
    }    
  }
  printf("\n\n \tRestrições:\n\n");
  for (int i = 0; i < qtyRests; i++ ){
    double *restriction = *(rawRestrictions + i);
    printf("\t");
    for (int j = 0; j< qtyTerms; j++ ){
      //DEBUG_PRINT(("\n i: %d j: %d  restrictionAddr: %p restrictionTermAddr/Value: %p %f String: ", i, j, (restriction), (restriction + j), *(restriction + j) ));
      if ( j < qtyTerms - 2 ){
        printf(" %.0fX%i", *(restriction + j), j+1 );
        printf(" + ");
      } else if(j == (qtyTerms - 2)) {
        printf(" %.0fX%i", *(restriction + j), j+1 );
      } else {
        printf(" <= %.0f\n", *(restriction + j));
      }   
    }
  }
  for (int i=1; i<=qtyVar; i++){
    printf("\t\tX%i >= 0\n", i);
  }

}

void populate_collumn_names(){

  int zbtermsSize=2;
  int totalSize = qtyVariables + qtySlackVariables + qtyArtificialVariables + zbtermsSize;
  int zLimit = 0;
  int xnLimit = qtyVariables;
  int slackLimit = qtyVariables + qtySlackVariables;
  int xArtLimit = totalSize-2;
  collumnNames = malloc(sizeof(char *) * (totalSize) );


  for (int i =0; i < totalSize; i++){
    char *collumnName = malloc(sizeof(char) * VARNAME_SIZE);
    if (i == zLimit) {
      char item[] = "Z";
      strcat(collumnName, item);     
    } else if (i <= xnLimit ){
      char prefix[]="X";
      strcat(collumnName, prefix);
      char item[VARNAME_SIZE-2];
      sprintf(item, "%d" ,i);
      strcat(collumnName, item);
    } else if (i <= slackLimit ) {
      strcat(collumnName, "X");
      char index[VARNAME_SIZE-2];
      sprintf(index, "%d" ,i);
      strcat(collumnName, index);
      strcat(collumnName, "f");
    } else if (i <= xArtLimit ){
      strcat(collumnName, "X");
      char item[VARNAME_SIZE-2];
      sprintf(item, "%d" ,i);
      strcat(collumnName, item);
      strcat(collumnName, "a");
    } else{
      strcat(collumnName, "b");  
    }
    *(collumnNames + i) = collumnName;
    //DEBUG_PRINT(("%s\t%i\t%p\t%p\n", collumnName, strlen(collumnName), collumnName, *(collumnNames + i)));
  
  }
}

double * calculate_pivot_line( int lineIndex, int CollIndex){

  DEBUG_PRINT(( "\n\n//// CALC PIVOT LINE //// lineIndex: %i, CollIndex %i\n\n", lineIndex, CollIndex ));

  double *currentLine = *(tableau + lineIndex);
  double pivot = *(currentLine + CollIndex);
  double *pivotLine = malloc(sizeof(double)* qtyCollumnsTableau);
  for (int i=0; i< qtyCollumnsTableau ;i++){
    double oldValue = *(currentLine + i);
    double newValue = oldValue / pivot;
    *(pivotLine + i) = newValue;

    DEBUG_PRINT(( "\n i: %i, oldValue: %f, pivot: %f, newValue: %f  \n", i, oldValue, pivot, newValue ));  
  }

  return pivotLine;
}


double * calculate_new_tableau_line( double *pivotLine, int lineIndex,  int CollIndex){

  DEBUG_PRINT(( "\n\n//// CALC NEW TABLEAU LINE //// lineIndex: %i, CollIndex %i\n\n", lineIndex, CollIndex ));

  double *currentLine = *(tableau + lineIndex);
  double *newLine = malloc(sizeof(double) * qtyCollumnsTableau);
  double invMultiplier = (*(currentLine + CollIndex) * -1.0);
  for (int i=0 ; i< qtyCollumnsTableau; i++){
    double pivotValue = *(pivotLine + i);
    double oldValue = *(currentLine + i);
    double newValue =  (pivotValue * invMultiplier) + oldValue;
    *(newLine + i) = newValue;
    DEBUG_PRINT(( "\n @@@@ Computed new Value: %f / Original: %f / pivotValue: %f / InvMultiplier: %f  @@@@ \n\n", newValue, oldValue, pivotValue, invMultiplier ));
  }
  return  newLine;
}

bool achieve_optimal_solution(){

    bool isOptimalSolution = true;
    double *zLine = *(tableau);

    for (int i=0; i< qtyCollumnsTableau; i++){
      if(*(zLine + i) < 0){
        isOptimalSolution = false;
      }
    } 

    return isOptimalSolution;
}

void update_tableau_line(double *line, int lineIndex){

  *(tableau + lineIndex) = line;

}

void update_line_names(int outterIndex, int innerIndex){

  DEBUG_PRINT(( "\n\n//// UPDATE LINE NAMES //// outterIndex: %i, InnerIndex %i\n\n", outterIndex, innerIndex ));

  if ((outterIndex >= 0 && outterIndex < qtyLinesTableau) && (innerIndex >=0 && innerIndex < qtyCollumnsTableau)){
      
      DEBUG_PRINT(( "\n @@@@ Replacing %s by %s @@@@ \n\n",  *(lineNames + outterIndex),  *(collumnNames + innerIndex) ));

      char *innerVarName = *(collumnNames + innerIndex);
      *(lineNames + outterIndex) = innerVarName;
  }else{
    DEBUG_PRINT(( "\n\t ###### ERROR: at least one INDEX is out of tableau bounds. \n\n",  *(lineNames + outterIndex),  *(collumnNames + innerIndex) ));
  }

}


void populate_line_names (){

  int qtyLines = qtySlackVariables;
  int qtyBasicVariables = qtyVariables;
  int xnfOffset = qtyBasicVariables;
  lineNames = malloc(sizeof(char *) * (qtyLines));
  for (int i = 0 ; i < qtyLinesTableau ; i++){
    char *lineItem = malloc(sizeof(char)*VARNAME_SIZE);
    if (i==0) { 
      strcat(lineItem, "Zline");
    } else{
      strcat(lineItem,"X");
      char index[VARNAME_SIZE-1];
      //DEBUG_PRINT("%i %i %i",i + xnfOffset, i , xnfOffset );
      int j = i + xnfOffset;
      sprintf(index, "%i", j);
      strcat(lineItem, index);
      strcat(lineItem, "f");
    }
    *(lineNames + i) = lineItem;
    //DEBUG_PRINT(("--\n%i%s\t%i\t%p\t%p\n", i, lineItem, strlen(lineItem), lineItem, *(lineNames + i)));
    char *lineItem2 = *(lineNames + i);
    //DEBUG_PRINT(("==\n%s\t%i\t%p\t%p\n", lineItem2, strlen(lineItem2), lineItem2, *(lineNames + i)));
  }

}

double ** populate_tableau(double *normalizedObjectiveFunction,double **normalizedRestrictions){

  int qtyLines = qtyLinesTableau;
  int qtyCollumns = qtyCollumnsTableau; 
  double **localTableau = malloc(sizeof(double *) * qtyLines) ;

  for (int i = 0; i< qtyLines; i++) {
    if (i==0) { 
      *(localTableau) = (normalizedObjectiveFunction);  
    } else {
      *(localTableau + i ) = *(normalizedRestrictions + (i-1));
    }
  }


  return localTableau; 

}

int detect_inner_var(){


  DEBUG_PRINT(( "\n\n//// DETECT INNER VAR //// \n\n"));

  int qtyTerms = qtyCollumnsTableau-1;
  int index = -1;
  double *zLine = *(tableau);
  
  double a = 0;
  for (int i=0; i < qtyTerms ; i++){
    double current = *(zLine + i);
    if( (current < 0) && (current < a) ){
      a = current;
      index = i;
      DEBUG_PRINT((" @@@@ Found new Coll: %i Value: %f @@@@", index, a));
    }
    DEBUG_PRINT (("\t\ta: %.0f, i: %d, index: %d, current: %.0f\n", a,i, index, current));
  }    

  DEBUG_PRINT((" >>>>>> Result Coll: %i Value of Coll: %f @@@@", index, a));

  return index;
}

int detect_outter_var( int innerIndex ){

  DEBUG_PRINT(( "\n\n//// DETECT OUTTER VAR //// innerIndex: %i\n\n", innerIndex ));

  int qtyLines = qtyLinesTableau;
  int bTermPos = qtyCollumnsTableau - 1;
  int outterRowIndex = -1;
  double lesserValue = FLT_MAX;
  for (int i = 0; i < qtyLines; i++ ){
    double *line = *(tableau + i);
    if ( i==0 ){
      // outterRowIndex = 0;
      // lesserValue = ;
    } else { 
      double innerValue =  *(line + i);
      double bValue = *(line + bTermPos);
      double computedValue = bValue / innerValue;
      if ( (computedValue < lesserValue) && (computedValue >= 0.0) ) {
        outterRowIndex = i;
        lesserValue = computedValue;
        DEBUG_PRINT((" @@@@ Found new row: %i Value: %f @@@@", outterRowIndex, lesserValue));
      }
      DEBUG_PRINT(( "\n i: %i, innerValue: %f, bValue: %f, computedValue: %f //// lesserValue: %f outterRowIndex %f \n", 
        i, innerValue, bValue, computedValue, lesserValue, outterRowIndex ));
    }
    
  }

  DEBUG_PRINT((" >>>>>> Result row: %i LesserValue: %f @@@@", outterRowIndex, lesserValue));
  return outterRowIndex;  
}


double * normalize_objective_function(double *rawObjectiveFunction){

  int qtyTerms = qtyCollumnsTableau;
  int qtyRawRestTerms = qtyVariables + 1; 


  double *normalizedObjectiveFunction = malloc(sizeof(double *) * qtyTerms);


  for (int i = 0; i < qtyTerms ; i++ ) {
    
    if (i == 0 ){
      *(normalizedObjectiveFunction + i ) = 1;
    } else if(i < qtyRawRestTerms){
      *(normalizedObjectiveFunction + i ) = *(rawObjectiveFunction + i) * (-1.0);
    } else {
      *(normalizedObjectiveFunction + i ) = 0;
    }
    //DEBUG_PRINT ((" item %i : %f  == %f \n", i , *(normalizedObjectiveFunction + i ), *(rawObjectiveFunction + i)  ));
  }


  return normalizedObjectiveFunction;

}

double ** normalize_restrictions(double **rawRestrictions){

  int qtyRules = qtyRestrictions;
  int qtyTerms = qtyCollumnsTableau;
  int qtyRawRestTerms = qtyVariables + 1;

  int bTermPos = qtyCollumnsTableau - 1;
  int rawBtermPos = qtyRawRestTerms - 1;
  
  

  double **normalizedRestrictions = malloc(sizeof(double *) * qtyRestrictions);

  for (int i=0; i < qtyRules; i++) {
    double *restriction = malloc (sizeof(double)*qtyTerms);
    double *rawRestriction = *(rawRestrictions + i);
    for (int j=0; j< (qtyTerms-1); j++){

      // se na regiao original copie, senao adicione 1 em posição da variavel de folga
      // adicione o rightvalue na ultima posição
      if (j == 0){
        *(restriction) = 0;  
      } else if  ( j == rawBtermPos){
        //copia termo b para o fim do array
        *(restriction + bTermPos) = *(rawRestriction + j);
        //preenche a posição corrente do array com o valor correspondente
        *(restriction + j ) = *(rawRestriction + (j-1));
      } else if( j < (rawBtermPos)){
        *(restriction + j ) = *(rawRestriction + (j-1));
        //DEBUG_PRINT (("bterm = %i, rawbterm = %i j = %i, new = %f , raw = %f\n", bTermPos, rawBtermPos, j, *(restriction + j), *(rawRestriction + (j-1))  ));
      } else if ( j == (qtyRawRestTerms + i) ){
        *(restriction + j) = 1;
      } else if ( j < bTermPos){
        *(restriction + j) = 0;        
      }
    }
    *(normalizedRestrictions + i) = restriction;
  }

  return normalizedRestrictions;
}

void solve(){

  print_tableau_equations();  
  print_tableau_state();

  do{
    iterations ++;
    int innerCollIndex = detect_inner_var();
    int outterRowIndex = detect_outter_var(innerCollIndex);
    double *pivotLine = calculate_pivot_line(outterRowIndex, innerCollIndex);
    update_line_names(outterRowIndex, innerCollIndex);
    for (int i=0; i<qtyLinesTableau; i++){
      if(i == outterRowIndex){
        update_tableau_line(pivotLine, outterRowIndex);
      }else{
        double *newLine = calculate_new_tableau_line(pivotLine, i, innerCollIndex);
        update_tableau_line(newLine,i);  
      }
    }


    print_tableau_state();
    if(max_iterations >=0){
      if((iterations >= max_iterations)){

        break;
      }
    }

  }while(!achieve_optimal_solution() );
  
  printf("\n\n########### Solution Achieved After  %i Iterations ###########\n\n",iterations);

  print_tableau_equations();
  print_tableau_state();
  print_tableau_variables();



}

/*
 * Inicializa o modelo enquanto o processamento de input não é habilitado
 *
 */
void setup_test_model_simplex_1(){

  qtyVariables = 2;

  qtyRestrictions = 2;
  qtySlackVariables = qtyRestrictions;
  qtyArtificialVariables = 0;

  qtyCollumnsTableau = 1 + qtyVariables + qtySlackVariables + qtyArtificialVariables + 1;
  qtyLinesTableau = 1 + qtyRestrictions;

  //Hand Crafted Values
  rawObjectiveFunction = malloc(sizeof(double *)*(qtyVariables +1));
  rawObjectiveFunction[0] = 0;
  rawObjectiveFunction[1] = 10;
  rawObjectiveFunction[2] = 12;




  rawRestrictions = malloc(sizeof(double *) * qtyRestrictions);
  double *rawRestriction1 = malloc( sizeof(double) * 3);
  *(rawRestriction1) = 1.0;
  *(rawRestriction1 + 1) = 1.0;
  *(rawRestriction1 + 2) = 100.0;

  double *rawRestriction2 = malloc (sizeof(double) * 3);
  *(rawRestriction2) = 1.0;
  *(rawRestriction2 + 1) = 3.0;
  *(rawRestriction2 + 2) = 270.0;

  //double rawRestriction3[3] = { 20, 30, 550};
  *(rawRestrictions) = rawRestriction1;
  *(rawRestrictions + 1) = rawRestriction2;
  //rawRestrictions[2] = rawRestriction3;

  populate_collumn_names();
  populate_line_names();


  normRestrictions = normalize_restrictions(rawRestrictions); 
  normObjectiveFunction = normalize_objective_function(rawObjectiveFunction);

  tableau = populate_tableau(normObjectiveFunction, normRestrictions);

}

/*
 * Inicializa o modelo enquanto o processamento de input não é habilitado
 *
 */
void setup_test_model_simplex_2(){

  qtyVariables = 2;

  qtyRestrictions = 4;
  qtySlackVariables = qtyRestrictions;
  qtyArtificialVariables = 0;

  qtyCollumnsTableau = 1 + qtyVariables + qtySlackVariables + qtyArtificialVariables + 1;
  qtyLinesTableau = 1 + qtyRestrictions;

  //Hand Crafted Values
  rawObjectiveFunction = malloc(sizeof(double *)*(qtyVariables +1));
  rawObjectiveFunction[0] = 0;
  rawObjectiveFunction[1] = 600;
  rawObjectiveFunction[2] = 800;




  rawRestrictions = malloc(sizeof(double *) * qtyRestrictions);
  double *rawRestriction1 = malloc( sizeof(double) * 3);
  *(rawRestriction1) = 1.0;
  *(rawRestriction1 + 1) = 1.0;
  *(rawRestriction1 + 2) = 100.0;

  double *rawRestriction2 = malloc (sizeof(double) * 3);
  *(rawRestriction2) = 3.0;
  *(rawRestriction2 + 1) = 2.0;
  *(rawRestriction2 + 2) = 240.0;

  double *rawRestriction3 = malloc (sizeof(double) * 3);
  *(rawRestriction3) = 1.0;
  *(rawRestriction3 + 1) = 0.0;
  *(rawRestriction3 + 2) = 60.0;

  double *rawRestriction4 = malloc (sizeof(double) * 3);
  *(rawRestriction4) = 0.0;
  *(rawRestriction4 + 1) = 1.0;
  *(rawRestriction4 + 2) = 80.0;

  //double rawRestriction3[3] = { 20, 30, 550};
  *(rawRestrictions) = rawRestriction1;
  *(rawRestrictions + 1) = rawRestriction2;
  *(rawRestrictions + 2) = rawRestriction3;
  *(rawRestrictions + 3) = rawRestriction4;
  //rawRestrictions[2] = rawRestriction3;

  populate_collumn_names();
  populate_line_names();


  normRestrictions = normalize_restrictions(rawRestrictions); 
  normObjectiveFunction = normalize_objective_function(rawObjectiveFunction);

  tableau = populate_tableau(normObjectiveFunction, normRestrictions);

}


/*
 * Inicializa o modelo enquanto o processamento de input não é habilitado
 *
 */
void setup_test_model_simplex_3(){

  max_iterations = 20;


  qtyVariables = 2;

  qtyRestrictions = 3;
  qtySlackVariables = qtyRestrictions;
  qtyArtificialVariables = 0;

  qtyCollumnsTableau = 1 + qtyVariables + qtySlackVariables + qtyArtificialVariables + 1;
  qtyLinesTableau = 1 + qtyRestrictions;

  //Hand Crafted Values
  rawObjectiveFunction = malloc(sizeof(double *)*(qtyVariables +1));
  rawObjectiveFunction[0] = 0;
  rawObjectiveFunction[1] = 2;
  rawObjectiveFunction[2] = 3;




  rawRestrictions = malloc(sizeof(double *) * qtyRestrictions);
  double *rawRestriction1 = malloc( sizeof(double) * 3);
  *(rawRestriction1) = 1.0;
  *(rawRestriction1 + 1) = 3.0;
  *(rawRestriction1 + 2) = 9.0;

  double *rawRestriction2 = malloc (sizeof(double) * 3);
  *(rawRestriction2) = -1.0;
  *(rawRestriction2 + 1) = 2.0;
  *(rawRestriction2 + 2) = 4.0;

  double *rawRestriction3 = malloc (sizeof(double) * 3);
  *(rawRestriction3) = 1.0;
  *(rawRestriction3 + 1) = 1.0;
  *(rawRestriction3 + 2) = 6.0;


  //double rawRestriction3[3] = { 20, 30, 550};
  *(rawRestrictions) = rawRestriction1;
  *(rawRestrictions + 1) = rawRestriction2;
  *(rawRestrictions + 2) = rawRestriction3;

  populate_collumn_names();
  populate_line_names();


  normRestrictions = normalize_restrictions(rawRestrictions); 
  normObjectiveFunction = normalize_objective_function(rawObjectiveFunction);

  tableau = populate_tableau(normObjectiveFunction, normRestrictions);

}


/*
 * Inicializa o modelo enquanto o processamento de input não é habilitado
 *
 */
void setup_test_model_simplex_4(){

  max_iterations = 20;

  qtyVariables = 2;

  qtyRestrictions = 2;
  qtySlackVariables = qtyRestrictions;
  qtyArtificialVariables = 0;

  qtyCollumnsTableau = 1 + qtyVariables + qtySlackVariables + qtyArtificialVariables + 1;
  qtyLinesTableau = 1 + qtyRestrictions;

  //Hand Crafted Values
  rawObjectiveFunction = malloc(sizeof(double *)*(qtyVariables +1));
  rawObjectiveFunction[0] = 0;
  rawObjectiveFunction[1] = 120;
  rawObjectiveFunction[2] = 60;




  rawRestrictions = malloc(sizeof(double *) * qtyRestrictions);
  double *rawRestriction1 = malloc( sizeof(double) * 3);
  *(rawRestriction1) = 1.0;
  *(rawRestriction1 + 1) = 1.5;
  *(rawRestriction1 + 2) = 50.0;

  double *rawRestriction2 = malloc (sizeof(double) * 3);
  *(rawRestriction2) = 4.0;
  *(rawRestriction2 + 1) = 3.0;
  *(rawRestriction2 + 2) = 360.0;


  *(rawRestrictions) = rawRestriction1;
  *(rawRestrictions + 1) = rawRestriction2;

  populate_collumn_names();
  populate_line_names();


  normRestrictions = normalize_restrictions(rawRestrictions); 
  normObjectiveFunction = normalize_objective_function(rawObjectiveFunction);

  tableau = populate_tableau(normObjectiveFunction, normRestrictions);

}


void setup_test_model_simplex_5(){

  max_iterations = 20;

  qtyVariables = 2;

  qtyRestrictions = 2;
  qtySlackVariables = qtyRestrictions;
  qtyArtificialVariables = 0;

  qtyCollumnsTableau = 1 + qtyVariables + qtySlackVariables + qtyArtificialVariables + 1;
  qtyLinesTableau = 1 + qtyRestrictions;

  //Hand Crafted Values
  rawObjectiveFunction = malloc(sizeof(double *)*(qtyVariables +1));
  rawObjectiveFunction[0] = 0;
  rawObjectiveFunction[1] = 120;
  rawObjectiveFunction[2] = 60;




  rawRestrictions = malloc(sizeof(double *) * qtyRestrictions);
  double *rawRestriction1 = malloc( sizeof(double) * 3);
  *(rawRestriction1) = 1.0;
  *(rawRestriction1 + 1) = 1.5;
  *(rawRestriction1 + 2) = 50.0;

  double *rawRestriction2 = malloc (sizeof(double) * 3);
  *(rawRestriction2) = 4.0;
  *(rawRestriction2 + 1) = 3.0;
  *(rawRestriction2 + 2) = 360.0;


  *(rawRestrictions) = rawRestriction1;
  *(rawRestrictions + 1) = rawRestriction2;

  populate_collumn_names();
  populate_line_names();


  normRestrictions = normalize_restrictions(rawRestrictions); 
  normObjectiveFunction = normalize_objective_function(rawObjectiveFunction);

  tableau = populate_tableau(normObjectiveFunction, normRestrictions);

}


int main(){

  DEBUG_PRINT(( "DEBUG Habilitado \n"));
  //init_ncurses();  

  printf("\n\n >>>>>>>>>>>>>>>>>>>>> Simplex Solver <<<<<<<<<<<<<<<<<<<<<\n\n");

  printf("\n\t Insira o número do problema que deseja resolver [1-5]: ");
  int *problemId = malloc(sizeof(int));
  scanf("%i", problemId );

  printf("\n Loading problem %i ... \n", *problemId);

  switch (*problemId){

    case 1:
      setup_test_model_simplex_1();
      break;
    case 2:
      setup_test_model_simplex_2();
      break;
    case 3:
      setup_test_model_simplex_3();
      break;
    case 4:
      setup_test_model_simplex_4();
      break;
    case 5:
      setup_test_model_simplex_5();
      break;
    default:
      printf("\n Invalid problem %i, Loading problem 3 instead ... \n", *problemId);
      setup_test_model_simplex_3();
      break;      
  }

  solve();

  return 0;  
} 
