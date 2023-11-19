/* This is a sample Advection solver in C 
The advection equation-> \partial q / \partial t - u \cdot \nabla q(x,y) = 0
The grid of NX by NX evenly spaced points are used for discretization.  
The first and last points in each direction are boundary points. 
Approximating the advection operator by 1st order finite difference. 
*/
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include "advection.h"
#define BUFSIZE 512
/* ************************************************************************** */
int main ( int argc, char *argv[] ){
  if(argc!=2){
    printf("Usage: ./levelSet input.dat\n");
    return -1;  
  }
  static int frame=0;

  // Create an advection solver
  solver_t advc; 
  // Create uniform rectangular (Cartesian) mesh
  advc.msh = createMesh(argv[1]); 
  // Create time stepper 
  tstep_t tstep = createTimeStepper(advc.msh.Nnodes); 
  // Create Initial Field
  initialCondition(&advc);

  // Read input file for time variables 
  tstep.tstart = readInputFile(argv[1], "TSART");
  tstep.tend   = readInputFile(argv[1], "TEND");
  tstep.dt     = readInputFile(argv[1], "DT");
  tstep.time = 0.0; 

  // adjust time step size 
  int Nsteps = ceil( (tstep.tend - tstep.tstart)/tstep.dt);
  tstep.dt = (tstep.tend - tstep.tstart)/Nsteps;
  // Read input file for OUTPUT FREQUENCY i.e. in every 1000 steps
  int Noutput = readInputFile(argv[1], "OUTPUT_FREQUENCY");


  // write the initial solution i.e. q at t = tstart
  {
    char fname[BUFSIZ];
    sprintf(fname, "test_%04d.csv", frame++);
    solverPlot(fname, &advc.msh, advc.q);
  }

    clock_t start, end;
    double cpu_time_used;
  // ********************Time integration***************************************/
  // for every steps
  start = clock();
  for(int step = 0; step<Nsteps; step++){

     
    
    for(int stage=0; stage<tstep.Nstage; stage++){
      // Call integration function

      RhsQ(&advc, &tstep, stage); 

    }

    tstep.time = tstep.time+tstep.dt;

    if(step%Noutput == 0){
      char fname[BUFSIZ];
      sprintf(fname, "test_%04d.csv", frame++);
      solverPlot(fname, &advc.msh, advc.q);
      end = clock();
      cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
      printf("%lf\n",cpu_time_used);
      cpu_time_used=cpu_time_used-cpu_time_used;
      start = clock();
    }

  }
}

/* ************************************************************************** */
void RhsQ(solver_t *solver, tstep_t *tstep, int stage){

mesh_t *msh = &solver->msh;
for(int j=0; j<msh->NY; j++){
    for(int i=0; i<msh->NX; i++){
      
    int n = j * msh->NX + i;
    double DuqDx;
    double DvqDy;
    if (solver->u[n]>=0){
    DuqDx=(solver->u[n]*solver->q[n]-solver->u[msh->N2N[4*n]]*solver->q[msh->N2N[4*n]])/(msh->x[n]-msh->x[msh->N2N[4*n]]);
    }else{ 
    DuqDx=(-solver->u[n]*solver->q[n]+solver->u[msh->N2N[4*n+2]]*solver->q[msh->N2N[4*n+2]])/(-msh->x[n]+msh->x[msh->N2N[4*n+2]]);
    }
    int v_n=n+msh->Nnodes;
    if (solver->u[v_n]>=0){
    DvqDy=(solver->u[v_n]*solver->q[v_n]-solver->u[msh->N2N[4*n+3]+msh->Nnodes]*solver->q[msh->N2N[4*n+3]+msh->Nnodes])/(msh->y[n]-msh->y[msh->N2N[4*n+3]]);
    }else{ 
    DvqDy=(-solver->u[v_n]*solver->q[v_n]+solver->u[msh->N2N[4*n+1]+msh->Nnodes]*solver->q[msh->N2N[4*n+1]+msh->Nnodes])/(-msh->y[n]+msh->y[msh->N2N[4*n+1]]);
    }
    double rhsq=-(DuqDx+DvqDy);
    tstep->resq[n] = tstep->rk4a[stage] * tstep->resq[n] +tstep->dt * rhsq;
    solver->q[n] = solver->q[n] + tstep->rk4b[stage] * tstep->resq[n];
    }
  }
}

/* ************************************************************************** */
void initialCondition(solver_t *solver){
  mesh_t *msh = &(solver->msh); 

  solver->q = (double *)malloc(msh->Nnodes*sizeof(double)); 
  solver->u = (double *)malloc(2*msh->Nnodes*sizeof(double));

  for(int j=0; j<msh->NY; j++){
    for(int i=0; i<msh->NX; i++){
      double x_c=0.5;
      double y_c=0.75;
      double r=0.15;
      int n=msh->NX*j+i;
      solver->u[n] = sin( ( 0.5 + msh->x[n] ) * 3.14 * 4 ) * sin( ( 0.5 + msh->y[n] ) * 3.14 * 4 );
      solver->u[n+msh->Nnodes] = cos( ( 0.5 + msh->y[n] ) * 3.14 * 4 ) * cos( ( 0.5 + msh->x[n] ) * 3.14 * 4 );
      solver->q[n] = sqrt( pow( ( msh->x[n] - x_c ) , 2 ) + pow( ( msh->y[n] - y_c ) , 2 ) ) - r;
    }
  }

}



/* ************************************************************************** */
// void createMesh(struct mesh *msh){
mesh_t createMesh(char* inputFile){

  mesh_t msh; 

  msh.NY   = readInputFile(inputFile, "NY");
  msh.NX   = readInputFile(inputFile, "NX");
  msh.xmin   = readInputFile(inputFile, "XMIN");
  msh.xmax   = readInputFile(inputFile, "XMAX");
  msh.ymin   = readInputFile(inputFile, "YMIN");
  msh.ymax   = readInputFile(inputFile, "YMAX");

  msh.Nnodes = msh.NX*msh.NY;
  msh.x = (double *) malloc(msh.Nnodes*sizeof(double));
  msh.y = (double *) malloc(msh.Nnodes*sizeof(double));

  double dx = ( msh.xmax - msh.xmin ) / ( msh.NX - 1 );
  double dy = ( msh.xmax - msh.xmin ) / ( msh.NX - 1 );
  for(int j=0; j<msh.NY; j++){
    for(int i=0; i<msh.NX; i++){
    int n = j*msh.NX+i;
    msh.x[n] = msh.xmin + i * dx;
    msh.y[n] = msh.ymin + j * dy;
    }
  }

  // Create connectivity and periodic connectivity
  /* 
  for every node 4 connections east north west and south
  Nothe that periodic connections require specific treatment
  */
  msh.N2N = (int *)malloc(4*msh.Nnodes*sizeof(int));

  for(int j=0; j<msh.NY; j++){
    for(int i=0; i<msh.NX; i++){
      int n = j * msh.NX + i;
      int m = i * msh.NY + j;
      msh.N2N[4*n]=   n - i + ( n - 1 + msh.NX ) % msh.NX ; // West
      msh.N2N[4*n+1]= m - j + ( j + 1 ) % msh.NY + ( msh.NX - 1 ) * ( ( j + 1 ) % msh.NY ) - ( msh.NY - 1 ) * i; // North
      msh.N2N[4*n+2]= n - i + ( i + 1 ) % msh.NX; // East
      msh.N2N[4*n+3]= m - j + ( m - 1 + msh.NY ) % msh.NY + ( msh.NX - 1 ) * ( ( j - 1 + msh.NY ) % msh.NY ) - ( msh.NY - 1 ) * i; // South
      
    }
  }

  return msh; 
}

/* ************************************************************************** */
void solverPlot(char *fileName, mesh_t *msh, double *Q){
    FILE *fp = fopen(fileName, "w");
    if (fp == NULL) {
        printf("Error opening file\n");
        return;
    }

    fprintf(fp, "X,Y,Z,Q \n");
    for(int n=0; n< msh->Nnodes; n++){
      fprintf(fp, "%.8f, %.8f,%.8f,%.8f\n", msh->x[n], msh->y[n], 0.0, Q[n]);
    } 
}

/* ************************************************************************** */
double readInputFile(char *fileName, char* tag){
  FILE *fp = fopen(fileName, "r");
  if (fp == NULL) {
    printf("Error opening the input file\n");
    return -1;
  }
  char name[100];
  double in_val[0];
    while (fgets(name, 50, fp) != NULL) {
    int nline_index = strcspn(name, "\n");
    int bracket_index = strcspn(name, "[");
    if (2 == strlen(name)) {
        name[nline_index] = '\0';
    }else if(bracket_index==0){
        int checker = 1;
        name[strlen(name)-2]='\0';
        for (int i = 0; name[i] != '\0'; ++i) {
          
          if (tag[i] == name[i + 1] && checker==1){
            checker=1;
          }else{
            checker=0;
          };

        }
        if (checker==1) {
            fscanf(fp,"%lf",in_val);
            fclose(fp);
            return in_val[0];
        }
        

    }else{
      continue;
    }
    }

    printf("Error finding input tag in the input file\n");
    fclose(fp);
    return -1;
}


/* ************************************************************************** */
// Time stepper clas RK(4-5)
// resq = rk4a(stage)* resq + dt*rhsq
//  q = q + rk4b(stage)*resq
tstep_t createTimeStepper(int Nnodes){
  tstep_t tstep; 
  tstep.Nstage = 5; 
  tstep.resq = (double *)calloc(Nnodes,sizeof(double)); 
  tstep.rhsq = (double *)calloc(Nnodes,sizeof(double));
  tstep.rk4a = (double *)malloc(tstep.Nstage*sizeof(double));
  tstep.rk4b = (double *)malloc(tstep.Nstage*sizeof(double));
  tstep.rk4c = (double *)malloc(tstep.Nstage*sizeof(double));

  tstep.rk4a[0] = 0.0; 
  tstep.rk4a[1] = -567301805773.0/1357537059087.0; 
  tstep.rk4a[2] = -2404267990393.0/2016746695238.0;
  tstep.rk4a[3] = -3550918686646.0/2091501179385.0;
  tstep.rk4a[4] = -1275806237668.0/842570457699.0;
        
  tstep.rk4b[0] = 1432997174477.0/9575080441755.0;
  tstep.rk4b[1] = 5161836677717.0/13612068292357.0; 
  tstep.rk4b[2] = 1720146321549.0/2090206949498.0;
  tstep.rk4b[3] = 3134564353537.0/4481467310338.0;
  tstep.rk4b[4] = 2277821191437.0/14882151754819.0;
             
  tstep.rk4c[0] = 0.0;
  tstep.rk4c[1] = 1432997174477.0/9575080441755.0;
  tstep.rk4c[2] = 2526269341429.0/6820363962896.0;
  tstep.rk4c[3] = 2006345519317.0/3224310063776.0;
  tstep.rk4c[4] = 2802321613138.0/2924317926251.0;
  return tstep; 
}


