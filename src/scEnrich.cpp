// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
#include <vector>
using namespace Rcpp;


// ## One draws balls *without* replacement. At the end, one sees the result:
//   ## X - number of white balls drawn: # ID targets hit
//   ## M - number of white balls in the urn from the beginning: # ID targets
//   ## N - number of black balls in the urn from the beginning: # non-ID targets
//   ## K - number of draws made: # all targets hit
//   X <- n.id.targets.hit
//      M <- n.id.targets
// 	N <- n.non.id.targets
// 	   K <- n.all.targets.hit

// hyper.p.val <- 1-phyper(X-1,M,N,K)
// p.val.f <- fisher.test(matrix(c(X,M-X,K-X,N-(K-X)),nrow=2),
// 		       alternative="greater")$p.value


int cellRTFpairingCpp(LogicalVector cellExp, CharacterVector gList,
		      CharacterVector idTFset, List pthwGenes){
  int n = pthwGenes.size();
  double pval;
  int X,M,N,K;
  CharacterVector curPthwGenes;
  for(int i = 0; i < n; i++){
    curPthwGenes = as<CharacterVector>(pthwGenes[i]);
    //pval = 1.0 - phyper(X-1,M,N,K);
  }
}

std::unordered_set<std::string>& makeSet(std::unordered_set<std::string>& Set,
					 CharacterVector cV){
  int n = cV.size();
  
  for(int i = 0; i < n; i++)
    Set.emplace(cV[i]);
  return Set;
}


std::unordered_set<std::string>::const_iterator
get_set_idx(std::unordered_set<std::string> S,
	    std::string gname){
  return(S.find(gname));
}


std::unordered_set<std::string>& set_from_cvector(CharacterVector genes,
						  std::unordered_set<std::string>& S){
  size_t n = genes.size();
  S.reserve(n);
  for(int i = 0;i < n; i++)
    S.emplace(genes[i]);
  return(S);
}

bool genes_in_set(CharacterVector genes, CharacterVector set){
  if(na_omit(match(genes,set)).size() > 0)
    return true;
  else
    return false;
}

CharacterVector c_genes_in_set(CharacterVector genes, CharacterVector set){
  return (na_omit(match(genes,set)).size());
}

bool any_in_set(std::unordered_set<std::string> Sfrom,
		     std::unordered_set<std::string> Sto){
  for(auto it = Sfrom.begin(); it != Sfrom.end(); it++)
    if(Sto.find(*it) != Sto.end())
      return(true);
  return(false);
}


// [[Rcpp::export]]
List scRTFpairingCpp2(List pthwGenes, CharacterVector pthwNames, int cellIdx,
		     LogicalVector expArray, CharacterVector gList_,
		     CharacterVector TFLIST_, CharacterVector RLIST_, DataFrame TFTF_){

  int EP,NeP,EnP,NnP;//pthw enrichment counts, e.g. Exp in Pthw(EP), Not exp Not in Pthw (NnP) etc.
  //Form sets from vectors
  std::unordered_set<std::string> RLIST;
  std::unordered_set<std::string> TFLIST;
  std::unordered_set<std::string> gList;
  CharacterVector  TFTFSource, TFTFTarget, TFTFEffect;
  RLIST = set_from_cvector(RLIST_,RLIST);
  TFLIST = set_from_cvector(TFLIST_,TFLIST);
  gList = set_from_cvector(gList_,gList);
  TFTFSource = TFTF_["Source"];
  TFTFTarget = TFTF_["Target"];
  TFTFEffect = TFTF_["Effect"];
  // std::cout << "TFTF:" << TFTFSource[0] << "," << TFTFSource[10] << std::endl;
  // std::cout << "TFTF:" << TFTFTarget[10] << "," << TFTFTarget[100] << std::endl;
  // std::cout << "TFTF:" << TFTFEffect[10] << "," << TFTFEffect[100] << std::endl;
  std::unordered_set<std::string> Rexpd, TFexpd;
  Function fisher_test("fisher.test");
  Function p_adjust("p.adjust");

  List R_list_ext, r_lst_pthw;

  // Find expressed TF's and R's
  for(int g = 0; g < expArray.size(); g++){
    if(expArray[g] == true){
      if(RLIST.count(as<std::string>(gList_[g])))
	Rexpd.emplace(as<std::string>(gList_[g]));
      if(TFLIST.count(as<std::string>(gList_[g])))
	TFexpd.emplace(as<std::string>(gList_[g]));
    }
  }
  // std::cout << "Size of Rexpd:" << Rexpd.size() << std::endl;
  // std::cout << "Size of TFexpd:" << TFexpd.size() << std::endl;
  // std::cout << "Size of TFLIST:" << TFLIST.size() << std::endl;
  // std::cout << "Size of RLIST:" << RLIST.size() << std::endl;

  // Iterate over the pathways
  for(int p = 0; p < pthwGenes.size(); p++){
    std::unordered_set<std::string> currPthwGenes;
    currPthwGenes = set_from_cvector(as<CharacterVector>(pthwGenes[p]),currPthwGenes);
    std::unordered_set<std::string> pthwExpdR,pthwExpdTF;
    EP = EnP = NeP = NnP = 0;
    // Proceed only if any R and TF is expressed in the pathway
    if(any_in_set(Rexpd,currPthwGenes) &&
       any_in_set(TFexpd,currPthwGenes)){
      // Iterate over genes to find enrichment
      for(int g = 0; g < expArray.size(); g++){
	// Fisher test preparations
	if(expArray[g] == true){
	  if(currPthwGenes.count(as<std::string>(gList_[g]))){
	    EP += 1; // expressed in the pathway
	    // Form pthwExpdR and pthwExpdTF at the same run
	    if(Rexpd.count(as<std::string>(gList_[g])))
	      pthwExpdR.emplace(as<std::string>(gList_[g]));
	    if(TFexpd.count(as<std::string>(gList_[g])))
	      pthwExpdTF.emplace(as<std::string>(gList_[g]));
	  }
	  else
	    EnP += 1; // expressed not in the pathway
	}
	else {
	  if(currPthwGenes.count(as<std::string>(gList_[g])))
	    NeP += 1; // not expressed in the pathway
	  else
	    NnP += 1; // not expressed not in the pathway
	}
      }
      // Accounting for pthw genes not available in data => considered not expressed
      for(auto it = currPthwGenes.begin(); it != currPthwGenes.end(); it++)
	if(gList.count(*it) == 0)
	  NeP += 1;
      //Fisher test
      IntegerVector v = {EP,NeP,EnP,NnP};
      v.attr("dim") = Dimension(2,2);
      // std::cout << "V:" << v << std::endl;
      List ft_res = fisher_test(v,_["alternative"]="greater");
      double pval = ft_res["p.value"];
      // std::cout << "Pval:" << pval << std::endl;
      if(pval < 0.05){ // not adjusted here!
	// std::cout << "pthwExpdR:";
	// for(auto it = pthwExpdR.begin(); it != pthwExpdR.end(); it++)
	//   std::cout << *it << ",";
	// std::cout << pthwExpdR.size() << std::endl;
	// std::cout << "pthwExpdTF:";
	// for(auto it = pthwExpdTF.begin(); it != pthwExpdTF.end(); it++)
	//   std::cout << *it << ",";
	// std::cout << pthwExpdTF.size() << std::endl;
	// Targets of the interface TF's
	std::vector<std::string> pthwExpdTF_v (pthwExpdTF.begin(), pthwExpdTF.end());
	//Number of compatible (effect=+1) and/or expressed targets
	std::vector<int> pthwExpdTF_numT(pthwExpdTF_v.size()), pthwExpdTF_numTexpd(pthwExpdTF_v.size());
	if(pthwExpdTF.size() == 0){
	  // std::cout << "No if.TF expressed, drop" << std::endl;
	  r_lst_pthw.push_back(NA_LOGICAL);
	  continue; // take the next pathway here
	}
	for(int v = 0; v < pthwExpdTF_v.size(); v++){
	  // For each pthwExpdTF (if.TF) find:
	  for(int i = 0; i < TFTFSource.size(); i++){
	    if(pthwExpdTF_v[v] == as<std::string>(TFTFSource[i])){// if if.TF is source
	      // 1. Effect is +1
	      if(TFTFEffect[i] == "1"){// if effect = +1
		// 2. Target is expressed
		for(int g = 0; g < expArray.size(); g++){//loop of expression array
		  if(expArray[g] == true){
		    if(as<std::string>(gList_[g]) == as<std::string>(TFTFTarget[i]) ){// if target is expressed
		      pthwExpdTF_numTexpd[v]++;
		      break;//do not look further in expression data
		    }
		  }
		}
		pthwExpdTF_numT[v]++; // no check for availability in data, just compatible in effect's sign
	      }
	    }
	  }
	}
	//Return the list
	r_lst_pthw.push_back(List::create(_["r"] = pthwExpdR, _["p"] = pval,
					  _["if.TF"] = pthwExpdTF_v, _["n.expd.t.TF"] = pthwExpdTF_numTexpd,
					  _["n.t.TF"] = pthwExpdTF_numT, _["cell.id"] = cellIdx));
      }
      else {// if pval >= 0.05
	r_lst_pthw.push_back(NA_LOGICAL);// NA_LOGICAL? Produces NA in R
      }
    }
    else { // if no expressed R and TF in a pathway
      r_lst_pthw.push_back(NA_LOGICAL);// NA_LOGICAL? Produces NA in R
    }
  } // End of the pthw loop
  // Assign pathway names
  r_lst_pthw.names() = pthwNames;

  return r_lst_pthw;
}


// [[Rcpp::export]]
List scRTFpairingCpp(List pthwGenes, CharacterVector pthwNames, CharacterVector cellNames,
		     LogicalMatrix expTable, CharacterVector idTFset_, CharacterVector gList_,
		     CharacterVector TFLIST_, CharacterVector RLIST_, DataFrame TFTF_){

  int EP,NeP,EnP,NnP;//pthw enrichment counts
  //Form sets from vectors
  std::unordered_set<std::string> RLIST;
  std::unordered_set<std::string> TFLIST;
  std::unordered_set<std::string> gList;
  std::unordered_set<std::string> idTFset;
  CharacterVector  TFTFSource, TFTFTarget, TFTFEffect;
  RLIST = set_from_cvector(RLIST_,RLIST);
  TFLIST = set_from_cvector(TFLIST_,TFLIST);
  gList = set_from_cvector(gList_,gList);
  idTFset = set_from_cvector(idTFset_,idTFset);
  TFTFSource = TFTF_["Source"];
  TFTFTarget = TFTF_["Target"];
  TFTFEffect = TFTF_["Effect"];
  // std::cout << "TFTF:" << TFTFSource[0] << "," << TFTFSource[10] << std::endl;
  // std::cout << "TFTF:" << TFTFTarget[10] << "," << TFTFTarget[100] << std::endl;
  // std::cout << "TFTF:" << TFTFEffect[10] << "," << TFTFEffect[100] << std::endl;
  std::unordered_set<std::string> Rexpd, TFexpd;
  Function fisher_test("fisher.test");
  Function p_adjust("p.adjust");

  List R_list_ext, r_lst_pthw;
  for(int c = 0; c < expTable.ncol(); c++){

    // Find expressed TF's and R's
    for(int g = 0; g < expTable.nrow(); g++){
      if(expTable(g,c) == true){
	if(RLIST.count(as<std::string>(gList_[g])))
	  Rexpd.emplace(as<std::string>(gList_[g]));
	if(TFLIST.count(as<std::string>(gList_[g])))
	  TFexpd.emplace(as<std::string>(gList_[g]));
      }
    }
    // std::cout << "Size of Rexpd:" << Rexpd.size() << std::endl;
    // std::cout << "Size of TFexpd:" << TFexpd.size() << std::endl;
    // std::cout << "Size of TFLIST:" << TFLIST.size() << std::endl;
    // std::cout << "Size of RLIST:" << RLIST.size() << std::endl;

    // Iterate over the pathways
    for(int p = 0; p < pthwGenes.size(); p++){
      std::unordered_set<std::string> currPthwGenes;
      currPthwGenes = set_from_cvector(as<CharacterVector>(pthwGenes[p]),currPthwGenes);
      std::unordered_set<std::string> pthwExpdR,pthwExpdTF;
      EP = EnP = NeP = NnP = 0;
      // Proceed only if any R and TF is expressed in the pathway
      if(any_in_set(Rexpd,currPthwGenes) &&
	 any_in_set(TFexpd,currPthwGenes)){
	// Iterate over genes to find enrichment
	for(int g = 0; g < expTable.nrow(); g++){
	  // Fisher test preparations
	  if(expTable(g,c) == true){
	    if(currPthwGenes.count(as<std::string>(gList_[g]))){
	      EP += 1; // expressed in the pathway
	      // Form pthwExpdR and pthwExpdTF at the same run
	      if(Rexpd.count(as<std::string>(gList_[g])))
		pthwExpdR.emplace(as<std::string>(gList_[g]));
	      if(TFexpd.count(as<std::string>(gList_[g])))
		pthwExpdTF.emplace(as<std::string>(gList_[g]));
	    }
	    else
	      EnP += 1; // expressed not in the pathway
	  }
	  else {
	    if(currPthwGenes.count(as<std::string>(gList_[g])))
	      NeP += 1; // not expressed in the pathway
	    else
	      NnP += 1; // not expressed not in the pathway
	  }
	}
	// Accounting for pthw genes not available in data => considered not expressed
	for(auto it = currPthwGenes.begin(); it != currPthwGenes.end(); it++)
	  if(gList.count(*it) == 0)
	    NeP += 1;
	//Fisher test
	IntegerVector v = {EP,NeP,EnP,NnP};
	v.attr("dim") = Dimension(2,2);
	// std::cout << "V:" << v << std::endl;
	List ft_res = fisher_test(v,_["alternative"]="greater");
	double pval = ft_res["p.value"];
	// std::cout << "Pval:" << pval << std::endl;
	if(pval < 0.05){ // not adjusted here!
	  // std::cout << "pthwExpdR:";
	  // for(auto it = pthwExpdR.begin(); it != pthwExpdR.end(); it++)
	  //   std::cout << *it << ",";
	  // std::cout << pthwExpdR.size() << std::endl;
	  // std::cout << "pthwExpdTF:";
	  // for(auto it = pthwExpdTF.begin(); it != pthwExpdTF.end(); it++)
	  //   std::cout << *it << ",";
	  // std::cout << pthwExpdTF.size() << std::endl;
	  CharacterVector targets, effects;
	  std::unordered_set<std::string> id_targets, non_id_targets,
	    expd_id_targets, expd_non_id_targets, all_up_targets, all_expd_targets;
	  // Targets of the interface TF's
	  for(auto it = pthwExpdTF.begin(); it != pthwExpdTF.end(); it++){
	    for(int i = 0; i < TFTFSource.size(); i++){
	      if(*it == as<std::string>(TFTFSource[i])){
	      	targets.push_back(TFTFTarget[i]);
	      	effects.push_back(TFTFEffect[i]);
	      }
	    }
	  }
	  // ID and non-ID targets
	  for(int g = 0; g < targets.size(); g++){
	    if(idTFset.count(as<std::string>(targets[g])) &
	       effects[g] == "1"){
	      id_targets.emplace(as<std::string>(targets[g]));
	    }
	    if(idTFset.count(as<std::string>(targets[g])) == 0 &
	       effects[g] == "1"){
	      non_id_targets.emplace(as<std::string>(targets[g]));
	    }
	    if(effects[g] == "1")// Auxillary variable
	      all_up_targets.emplace(as<std::string>(targets[g]));
	  }
	  // id_targets = unique(id_targets);
	  // non_id_targets = unique(non_id_targets);
	  // std::cout << "ID targets:" << id_targets.size() << std::endl;
	  // std::cout << "Non-ID targets:" << non_id_targets.size() << std::endl;
	  // Expressed ID and non-ID targets
	  for(int g = 0; g < expTable.nrow(); g++){
	    if(expTable(g,c) == true){
	      if(id_targets.count(as<std::string>(gList_[g])))
		expd_id_targets.emplace(as<std::string>(gList_[g]));
	      if(non_id_targets.count(as<std::string>(gList_[g])))
		expd_non_id_targets.emplace(as<std::string>(gList_[g]));
	      if(all_up_targets.count(as<std::string>(gList_[g])))
		all_expd_targets.emplace(as<std::string>(gList_[g]));
	    }
	  }
	  // std::cout << "Expd ID targets:" << expd_id_targets.size() << std::endl;
	  // std::cout << "Expd Non-ID targets:" << expd_non_id_targets.size() << std::endl;

	  // ========================
	  // Counts for measures
	  // ========================
	  // Vector equivalents of unordered sets, needed for sorting and set operations
	  std::vector<std::string> pthwExpdTF_v (pthwExpdTF.begin(), pthwExpdTF.end());
	  std::vector<std::string> idTFset_v (idTFset.begin(), idTFset.end());
	  std::vector<std::string> expd_id_targets_v (expd_id_targets.begin(), expd_id_targets.end());
	  std::vector<std::string> expd_non_id_targets_v (expd_non_id_targets.begin(), expd_non_id_targets.end());
	  std::vector<std::string> all_expd_targets_v(all_expd_targets.begin(), all_expd_targets.end());
	  std::vector<std::string> tmp_v, tmp_v2,idTFhit;
	  std::vector<std::string>::iterator it;
	  // Sort for the purposes of the set operations
	  std::sort(pthwExpdTF_v.begin(),pthwExpdTF_v.end());
	  std::sort(idTFset_v.begin(),idTFset_v.end());
	  std::sort(expd_id_targets_v.begin(),expd_id_targets_v.end());
	  std::sort(expd_non_id_targets_v.begin(),expd_non_id_targets_v.end());
	  std::sort(all_expd_targets_v.begin(),all_expd_targets_v.end());

	  // ========================
	  // Number of ID targets hit
	  // n_id_targets_hit = union_(expd_id_targets,intersect(pthwExpdTF,idTFset_)).size();
	  
	  //	  std::vector<std::string> idTFhit(pthwExpdTF.size()+idTFset.size());

	  // *** Intersection: Common genes in if.TF and id.TF sets
	  // prepare for the result
	  tmp_v.resize(pthwExpdTF.size() + idTFset.size());
	  // compute the intersection
	  it = std::set_intersection(pthwExpdTF_v.begin(),pthwExpdTF_v.end(),
				     idTFset_v.begin(),idTFset_v.end(),tmp_v.begin());
	  tmp_v.resize(it - tmp_v.begin()); // the result is sorted!
	  // *** Union: all genes belonging to the previous result and expd.id.targets
	  tmp_v2.resize(expd_id_targets.size() + tmp_v.size());
	  it = std::set_union(expd_id_targets_v.begin(),expd_id_targets_v.end(),
	  		      tmp_v.begin(),tmp_v.end(),tmp_v2.begin());
	  tmp_v2.resize(it - tmp_v2.begin());
	  int n_id_targets_hit = tmp_v2.size();
	  idTFhit = tmp_v2; // Keep the result for report purposes
	  // ============================
	  // Number of non-ID targets hit
	  // n_non_id_targets_hit = union_(expd_non_id_targets,setdiff(pthwExpdTF,idTFset_)).size();
	  
	  // *** Difference: select genes that differ among the two: pthwExpdTF and idTFset
	  // prepare the result
	  tmp_v.resize(pthwExpdTF.size()+idTFset.size());
	  // compute the difference
	  it = std::set_difference(pthwExpdTF_v.begin(),pthwExpdTF_v.end(),
	  			   idTFset_v.begin(),idTFset_v.end(),tmp_v.begin());
	  tmp_v.resize(it - tmp_v.begin());// the result is sorted!
	  // *** Union: all genes belonging to the previous result and expd.non.id.targets
	  tmp_v2.resize(expd_non_id_targets.size() + tmp_v.size());
	  it = std::set_union(expd_non_id_targets_v.begin(),expd_non_id_targets_v.end(),
	  		      tmp_v.begin(),tmp_v.end(),tmp_v2.begin());
	  tmp_v2.resize(it - tmp_v2.begin());
	  int n_non_id_targets_hit = tmp_v2.size();

	  // =============================
	  // Number of all targets hit
	  // n_all_targets_hit = union_(all_expd_targets,pthwExpdTF).size();
	  tmp_v.resize(all_expd_targets.size()+pthwExpdTF.size());
	  it = std::set_union(all_expd_targets_v.begin(),all_expd_targets_v.end(),
	  		      pthwExpdTF_v.begin(),pthwExpdTF_v.end(),tmp_v.begin());
	  tmp_v.resize(it - tmp_v.begin());
	  int n_all_targets_hit = tmp_v.size();

	  // std::cout << "Counts:" << n_id_targets_hit << "," << n_non_id_targets_hit <<
	  //   "," << n_all_targets_hit << std::endl;
	  if(n_all_targets_hit != (n_id_targets_hit + n_non_id_targets_hit))
	    std::cerr << "Warning: all targets hit do not match:" << n_all_targets_hit << "," <<
	      n_id_targets_hit << "," << n_non_id_targets_hit << ".\n";
	  int X = n_id_targets_hit;
	  int M = idTFset.size();
	  int N = TFLIST.size() - M;
	  int K = n_all_targets_hit;
	  int TP = n_id_targets_hit;
	  int FP = n_non_id_targets_hit;
	  double precision = TP/(1.0*(TP+FP));
	  double recall = TP/(1.0*M);
	  double F1 = 2.0*precision*recall/(1.0*(precision + recall));
	  double hyper_p_val = 1.0 - R::phyper(X-1,M,N,K,true,false);
	  // std::cout << "Stat:" << X << "," << M << "," << N << "," << K << "\n";
	  // std::cout << "Measures:" << TP << "," << FP << "," << precision << "," <<
	  //   recall << "," << F1 << "," << hyper_p_val << "\n";
	  //Return the list
	  CharacterVector idTFset_out(idTFhit.size());
	  for(int i = 0; i < idTFhit.size(); i++)
	    idTFset_out[i] = idTFhit[i];
	  r_lst_pthw.push_back(List::create(_["r"] = pthwExpdR, _["if.TF"] = pthwExpdTF,
	  				    _["id.TF"] = idTFset_out, _["p"] = pval,
					    _["hp"] = hyper_p_val, _["pr"] = precision,
					    _["rc"] = recall, _["f1"] = F1));
	}
	else {// if pval >= 0.05
	  r_lst_pthw.push_back(NA_LOGICAL);// NA_LOGICAL? Produces NA in R
	}
      }
      else { // if no expressed R and TF in a pathway
	r_lst_pthw.push_back(NA_LOGICAL);// NA_LOGICAL? Produces NA in R
      }
    }
    // Assign pathway names
    r_lst_pthw.names() = pthwNames;
    // Clear the vectors before the next cell is processed
    Rexpd.erase(Rexpd.begin(),Rexpd.end());
    TFexpd.erase(TFexpd.begin(),TFexpd.end());
    R_list_ext.push_back(r_lst_pthw);
    r_lst_pthw.erase(0,r_lst_pthw.size());
  }
  // Assign cell names
  R_list_ext.names() = cellNames;
  return R_list_ext;
}

