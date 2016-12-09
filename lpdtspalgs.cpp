/*******************************************************************************
 * MC658 - Projeto e Análise de Algoritmos III - 2s2016
 * Prof.: Flavio Keidi Miyazawa
 * PED: Mauro Henrique Mulati
 ******************************************************************************/

/***************************
 * Wendrey Lustosa Cardoso *
 * RA: 148234			   *
 ***************************/

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "lpdtspalgs.h"
#include "time.h"
#include "gurobi_c++.h"

bool naive(const LpdTspInstance &l, LpdTspSolution  &s, int tl);

//------------------------------------------------------------------------------

// Heurística construtiva para obter solução 

bool constrHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl) {

	clock_t st = clock();
	int si, ti;
	double cost, load = 0;
	bool done = false;
	DNode v, node;
	DNodeBoolMap inTour(l.g);
	vector <bool> carrying (l.k, false);
	LpdTspSolution sol;

	// Marca todos os vértices como não pertencentes ao tour
	
	for (DNodeIt n(l.g); n != INVALID; ++n) 
		inTour[n] = false;
		
	// Adiciona o vértice inicial ao tour

	inTour[l.depot] = true;
	sol.cost = 0.0;
	sol.tour.push_back(l.depot);
	
	// Adiciona os outros vértices enquanto o tour não estiver completo
	
	while (!done) {
	
		// Verifica restrição de tempo
		
		if (tl < (clock() - st) / CLOCKS_PER_SEC)
			return false;
		
		cost = -1;
	
		// Encontra o vizinho de menor custo para o qual seguir
	
		for (OutArcIt e(l.g, sol.tour.back()); e != INVALID; ++e) {

			v = l.g.target(e);
			
			// Verifica se o nó já está no tour

			if (inTour[v])
				continue;
			
			// Se precisar coletar um item em v, verifica se pode carregar o item
			// Se precisar entregar um item em v, verifica se possui o item
			// Para todos os vizinhos que pode seguir, guarda o de menor custo

			si = l.s[v];
			ti = l.t[v];

			if (si != 0) {
				if (l.items[si-1].w + load <= l.capacity) {
					if (l.weight[e] < cost || cost == -1) {
						node = v;
						cost = l.weight[e];
					}
				}
			}

			else if (ti != 0) {
				if (carrying[ti-1]) {
					if (l.weight[e] < cost || cost == -1) {
						node = v;
						cost = l.weight[e];
					}			
				}
			}
			
		}
		
		// Verifica se chegou a uma solução inviável
		
		if (cost == -1)
			return false;
			
		// Adiciona o vértice ao tour
	
		inTour[node] = true;		
		sol.tour.push_back(node);
		sol.cost += cost;

		// Coleta ou entrega o item referente ao vértice
			
		si = l.s[node];
		ti = l.t[node];	
			
		if (si != 0) {
			load += l.items[si-1].w;
			carrying[si-1] = true; 
		}
		
		if (ti != 0) {
			load -= l.items[ti-1].w;
			carrying[ti-1] = false;	
		}
			
		// Verifica se o tour está completo
			
		if (sol.tour.size() == (unsigned) l.n)
			done = true;	
			
	}
	
	// Verifica se é possível voltar ao depósito
	// Retorna a solução viável encontrada
	
	for (OutArcIt e(l.g, sol.tour.back()); e != INVALID; ++e) {
		if (l.g.target(e) == l.depot) {	
			s.cost = sol.cost + l.weight[e];
			s.tour = sol.tour;
		}
	}

	return false;

}

//------------------------------------------------------------------------------
bool metaHeur(const LpdTspInstance &l, LpdTspSolution  &s, int tl)
/* Implemente esta função, entretanto, não altere sua assinatura */
{
   return naive(l, s, tl);
}
//------------------------------------------------------------------------------
bool exact(const LpdTspInstance &l, LpdTspSolution  &s, int tl) {

	int i, j, k;
	double M = DBL_MAX;
	double lowerBound, upperBound;
	LpdTspSolution sol;

	// Associa um vertice a uma posicao

	k = 0;
	DNodeIntMap nodes(l.g);

	for (DNodeIt n(l.g); n != INVALID; ++n)
		nodes[n] = k++;		

	// Acha uma solução inicial com heurística construtiva

	constrHeur(l,sol,tl);

	// Inicializa o modelo

	lowerBound = 1;
	upperBound = sol.cost + 1;
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	model.set(GRB_StringAttr_ModelName, "LpdTsp");
	model.getEnv().set(GRB_DoubleParam_TimeLimit, tl);
	model.getEnv().set(GRB_DoubleParam_Cutoff, upperBound);
		
	// Ci é o custo das arestas para ir do depósito até o vértice i

	GRBVar* C = new GRBVar[l.n + 1];

	for (i = 0; i <= l.n; i++)
		C[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, "");

	// Aij é o peso dos itens carregados na aresta (i,j)

	GRBVar** A = new GRBVar*[l.n + 1];

	for (i = 0; i <= l.n; i++)
		A[i] = new GRBVar[l.n + 1];

	for (ArcIt e(l.g); e != INVALID; ++e)
		if (l.g.target(e) != l.depot)
			A[nodes[l.g.source(e)]][nodes[l.g.target(e)]] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_BINARY, "");
		else
			A[nodes[l.g.source(e)]][l.n] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_BINARY, "");

	// Xij = 1 se a aresta (i,j) é usada, Xij = 0 caso contrário

	GRBVar** X = new GRBVar*[l.n + 1];	
	
	for (i = 0; i <= l.n; i++)
		X[i] = new GRBVar[l.n + 1];
		
	for (ArcIt e(l.g); e != INVALID; ++e)
		if (l.g.target(e) != l.depot) 
			X[nodes[l.g.source(e)]][nodes[l.g.target(e)]] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_BINARY, "");
		else
			X[nodes[l.g.source(e)]][l.n] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_BINARY, "");			
		
	// Ui é auxiliar usada para que haja apena um tour
	
	GRBVar* U = new GRBVar[l.n + 1];
	
	for (i = 0; i <= l.n; i++)
		U[i] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_INTEGER, "");
		
	// (1): Sum de i = 1 até n de Aij - Sum de k = 1 até n de Ajk = Bj, para 1 <= j <= n 
	
	for (DNodeIt n(l.g); n != INVALID; ++n) {
		for (InArcIt in(l.g, n); in != INVALID; ++in) {
			for (OutArcIt out(l.g, n); out != INVALID; ++out) {
				i = nodes[l.g.source(in)];
				j = nodes[l.g.target(out)];
				k = nodes[n];
				GRBLinExpr expr = (2 - X[i][k] - X[k][j]) * M;
				if (l.s[n] > 0) {
					double b = l.items[l.s[n]-1].w;
					model.addConstr(A[k][j] - A[i][k] + expr >= b, "");
					model.addConstr(A[k][j] - A[i][k] <= b + expr, "");
				}
				else if (l.t[n] > 0) {
					double b = -l.items[l.t[n]-1].w;
					model.addConstr(A[k][j] - A[i][k] + expr >= b, "");
					model.addConstr(A[k][j] - A[i][k] <= b + expr, "");
				}
			}
		}
	}
	
	// (2): Aij <= Capacidade, para 1 <= i,j <= n 
	
	for (ArcIt e(l.g); e != INVALID; ++e) {
		i = nodes[l.g.source(e)];
		j = nodes[l.g.target(e)];
		model.addConstr(A[i][j] <= l.capacity * X[i][j], "");
	}
	
	for (InArcIt e(l.g, l.depot); e != INVALID; ++e)
		model.addConstr(A[nodes[l.g.source(e)]][l.n] == 0, "");

	for (OutArcIt e(l.g, l.depot); e != INVALID; ++e)
		model.addConstr(A[nodes[l.depot]][nodes[l.g.target(e)]] == 0, "");
	
	// (3): Cj >= Ci + Wij, para 1 <= i,j <= n
	
	for (ArcIt e(l.g); e != INVALID; ++e) {
		if (l.g.target(e) != l.depot) {
			GRBLinExpr expr = (1 - X[nodes[l.g.source(e)]][nodes[l.g.target(e)]]) * M;
			model.addConstr(C[nodes[l.g.target(e)]] + expr >= C[nodes[l.g.source(e)]] + l.weight[e], "");	
		}
		else {
			GRBLinExpr expr = (1 - X[nodes[l.g.source(e)]][l.n]) * M;
			model.addConstr(C[l.n] + expr >= C[nodes[l.g.source(e)]] + l.weight[e], "");				
		}
	}
			
	// (4): Ct > Cs, para todo par (s,t) de um item
	
	for (k = 1; k <= l.k; k++)
		model.addConstr(C[nodes[l.items[k].s]] <= C[nodes[l.items[k].t]], "");
	
	// (5): Sum de j = 1 até n de Xij = 1, para 1 <= i <= n

	for (DNodeIt n(l.g); n != INVALID; ++n) {
		GRBLinExpr expr = 0;
		for (OutArcIt out(l.g, n); out != INVALID; ++out)
			expr += X[nodes[n]][nodes[l.g.target(out)]];
		model.addConstr(expr == 1, "");
	}		
		
	// (6): Sum de i = 1 até n de Xij = 1, para 1 <= j <= n

	for (DNodeIt n(l.g); n != INVALID; ++n) {
		GRBLinExpr expr = 0;
		for (InArcIt in(l.g, n); in != INVALID; ++in) {
			if (n != l.depot) 
				expr += X[nodes[l.g.source(in)]][nodes[n]];
			else
				expr += X[nodes[l.g.source(in)]][l.n];		
		}
		model.addConstr(expr == 1, "");
	}		
	
	// (7): Ui - Uj + nXij <= n - 1, para i != j, 2 <= i,j <= n 
	
	for (i = 1; i <= l.n; i++)
		for (j = 1; j <= l.n; j++)
			if (i != j)
				model.addConstr(U[i] - U[j] + (l.n + 1) * X[i][j] <= l.n, "");
	
	// Objetivo: Min C
		
	GRBLinExpr obj = C[l.n];
	model.setObjective(obj, GRB_MINIMIZE);
	model.update();
	model.optimize();

	// Atribui solução

	s.tour.push_back(l.depot);
	
	while (s.tour.size() < (unsigned) l.n) 
		for (OutArcIt e(l.g, s.tour.back()); e != INVALID; ++e) 
			if (X[nodes[s.tour.back()]][nodes[l.g.target(e)]].get(GRB_DoubleAttr_X))
				s.tour.push_back(l.g.target(e));

	sol.lowerBound = model.get(GRB_DoubleAttr_ObjBound);
	sol.upperBound = model.get(GRB_DoubleAttr_ObjVal);
	sol.cost = upperBound;

	if (model.get(GRB_IntAttr_Status) == GRB_OPTIMAL)
		return true;	

	return false;

}

//------------------------------------------------------------------------------

bool naive(const LpdTspInstance &instance, LpdTspSolution  &sol, int tl)
/*
 * Algoritmo ingênuo para o LPD-TSP. Ideia:
 * constrNaiveHeur(l, s)
 *    s.tour.push_back(l.depot)
 *    while(s.tour.size() < 2*l.k+1)
 *       v = argmin_{v' in V} {d_{(v,v')} | (v' é adj a v) e ((v' é s) ou (v' é t de i cujo s é u em l.tour))}
 *       l.tour.push_back(v)
 */
{
   DNode v,
         vl;

   double vval,
          vlval;

   int i;

   sol.tour.clear();
   sol.cost = 0.0;

   v = instance.depot;
   sol.tour.push_back(v);

   while((int)sol.tour.size() < 2 * instance.k + 1 && v != INVALID){
      v    = INVALID;
      vval = DBL_MAX;

      for(OutArcIt o(instance.g, sol.tour.back()); o != INVALID; ++o){
         vl    = instance.g.target(o);
         vlval = DBL_MAX;

         i = 0;
         while(i < (int)sol.tour.size() && vl != sol.tour[i]) i++;
         if(i < (int)sol.tour.size()) continue;

         if(instance.s[vl] > 0){  // If DNode vl is start of an item
            vlval = instance.weight[o];
         }
         else if(instance.t[vl] > 0){  // If DNode vl is término of an item
            i = 0;
            while(i < (int)sol.tour.size() && instance.t[ vl ] != instance.s[ sol.tour[i] ]){  // Look for the start DNode of the item which terminates in DNode vl
               i++;
            }
            if(i < (int)sol.tour.size()){
               vlval = instance.weight[o];
            }
         }
         
         if(vlval < vval){
            v    = vl;
            vval = vlval;
         }
      }

      if(v != INVALID){
         sol.tour.push_back(v);
         sol.cost += vval;
      }
   }

   if(v == INVALID){
      sol.cost = DBL_MAX;
   }
   else{
      OutArcIt o(instance.g, sol.tour.back());
      for(; o != INVALID; ++o){
         if(instance.g.target(o) == sol.tour.front()) break;
      }
      if(o != INVALID){
         sol.cost += instance.weight[o];
      }
   }

   return false;
}
//------------------------------------------------------------------------------

