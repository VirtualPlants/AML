/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       StructureAnalysis: Identifying patterns in plant architecture and development
 *
 *       Copyright 1995-2019 CIRAD AGAP
 *
 *       File author(s): Yann Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id: stat_funs3.cpp 18820 2016-09-09 12:18:22Z guedon $
 *
 *       Forum for StructureAnalysis developers:
 *
 *  ----------------------------------------------------------------------------
 *
 *                      GNU General Public Licence
 *
 *       This program is free software; you can redistribute it and/or
 *       modify it under the terms of the GNU General Public License as
 *       published by the Free Software Foundation; either version 2 of
 *       the License, or (at your option) any later version.
 *
 *       This program is distributed in the hope that it will be useful,
 *       but WITHOUT ANY WARRANTY; without even the implied warranty of
 *       MERCHANTABILITY or FITNESS For A PARTICULAR PURPOSE. See the
 *       GNU General Public License for more details.
 *
 *       You should have received a copy of the GNU General Public
 *       License along with this program; see the file COPYING. If not,
 *       write to the Free Software Foundation, Inc., 59
 *       Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 *
 *  ----------------------------------------------------------------------------
 */



#include <iostream>
#include <string>
#include <vector>

#include "stat_tool/distribution.h"
#include "stat_tool/compound.h"
#include "stat_tool/convolution.h"
#include "stat_tool/discrete_mixture.h"
#include "stat_tool/distance_matrix.h"
#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/mixture.h"

#include "sequence_analysis/sequences.h"
#include "sequence_analysis/renewal.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"

#include "aml/ammodel.h"
#include "aml/parseraml.h"
#include "aml/kernel_err.h"
#include "aml/amstring.h"
#include "aml/array.h"

#include "aml/stat_module.h"


using namespace stat_tool;
using namespace sequence_analysis;


#define ERR_MSG_ARRAY STAT_err_msgs_aml
#define MODULE_NAME   "Statistics"


extern const char *STAT_err_msgs_aml[];

extern int nb_required_computation(const AMObjVector &args);
extern int* buildIntArray(const AMObjVector &args , int arg_index , const char *function ,
                          int output_index , int &nb_element);
extern double* buildRealArray(const AMObjVector &args , int arg_index , const char *function ,
                              int output_index , int &nb_element , bool filter);



/*--------------------------------------------------------------*
 *
 *  Conversion d'une loi empirique en loi.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToDistribution(const AMObjVector &args)

{
  DiscreteParametricModel *dist;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToDistribution"));

  CHECKCONDVA(args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToDistribution" ,
                          args[0].tag.string().data() , "FREQUENCY_DISTRIBUTION"));

  dist = ((DiscreteDistributionData*)((STAT_model*)args[0].val.p)->pt)->extract_model(error);

  if (dist) {
    STAT_model* model = new STAT_model(dist);
    return AMObj(AMObjType::DISTRIBUTION , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ToDistribution");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Conversion d'une loi en loi empirique.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToFrequencyDistribution(const AMObjVector &args)

{
  DiscreteDistributionData *histo;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToFrequencyDistribution"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTRIBUTION ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToFrequencyDistribution" ,
                          args[0].tag.string().data() , "DISTRIBUTION"));

  histo = ((DiscreteParametricModel*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

  if (histo) {
    STAT_model* model = new STAT_model(histo);
    return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ToFrequencyDistribution");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Pour compatibilite ascendante : Conversion d'une loi en loi empirique.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToHistogram(const AMObjVector &args)

{
  return STAT_ToFrequencyDistribution(args);
}

/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi elementaire,
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractDistribution(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractDistribution"));

  if (args[0].tag() == AMObjType::DISCRETE_MIXTURE) {
    string *pstr;
    DiscreteParametricModel *dist;
    const DiscreteMixture *mixture;
    const DiscreteMixtureData *mixture_data;
    StatError error;


    mixture = (DiscreteMixture*)((STAT_model*)args[0].val.p)->pt;
    mixture_data = mixture->get_mixture_data();

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Weight") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

      dist = new DiscreteParametricModel(*(mixture->get_weight()) ,
                                         (mixture_data ? mixture_data->get_weight() : 0));

      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }

    if (*pstr == "Component") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 3 ,
                              args[2].tag.string().data() , "INT"));

      dist = mixture->extract(error , args[2].val.i);

      if (dist) {
        STAT_model* model = new STAT_model(dist);
        return AMObj(AMObjType::DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Mixture") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

      dist = new DiscreteParametricModel(*((Distribution*)mixture) , (FrequencyDistribution*)mixture_data);

      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                "Weight or Component or Mixture");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::CONVOLUTION) {
    string *pstr;
    DiscreteParametricModel *dist;
    const Convolution *convolution;
    const ConvolutionData *convolution_data;
    StatError error;


    convolution = (Convolution*)((STAT_model*)args[0].val.p)->pt;
    convolution_data = convolution->get_convolution_data();

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Elementary") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 3 ,
                              args[2].tag.string().data() , "INT"));

      dist = convolution->extract(error , args[2].val.i);

      if (dist) {
        STAT_model* model = new STAT_model(dist);
        return AMObj(AMObjType::DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Convolution") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

      dist = new DiscreteParametricModel(*((Distribution*)convolution) , (FrequencyDistribution*)convolution_data);

      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                "Elementary or Convolution");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::COMPOUND) {
    string *pstr;
    DiscreteParametricModel *dist;
    const Compound *compound;
    const CompoundData *compound_data;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

    compound = (Compound*)((STAT_model*)args[0].val.p)->pt;
    compound_data = compound->get_compound_data();

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Sum") {
      dist = new DiscreteParametricModel(*(compound->get_sum_distribution()) ,
                                         (compound_data ? compound_data->get_sum_frequency_distribution() : 0));
    }
    else if (*pstr == "Elementary") {
      dist = new DiscreteParametricModel(*(compound->get_distribution()) ,
                                         (compound_data ? compound_data->get_frequency_distribution() : 0));
    }
    else if (*pstr == "Compound") {
      dist = new DiscreteParametricModel(*((Distribution*)compound) , (FrequencyDistribution*)compound_data);
    }
    else {
      genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                  "Sum or Elementary or Compound");
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(dist);
    return AMObj(AMObjType::DISTRIBUTION , model);
  }

  if (args[0].tag() == AMObjType::RENEWAL) {
    string *pstr;
    int time;
    renewal_distribution dist_type;
    DiscreteParametricModel *dist;
    StatError error;


    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "NbEvent") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 3));

      dist_type = NB_EVENT;

      if (args[2].tag() != AMObjType::INTEGER) {
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 3 ,
                    args[2].tag.string().data() , "INT");
        return AMObj(AMObjType::ERROR);
      }
      else {
        time = args[2].val.i;
      }
    }

    else {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

      if (*pstr == "InterEvent") {
        dist_type = INTER_EVENT;
      }
      else if (*pstr == "LengthBias") {
        dist_type = LENGTH_BIAS;
      }
      else if (*pstr == "Backward") {
        dist_type = BACKWARD_RECURRENCE_TIME;
      }
      else if (*pstr == "Forward") {
        dist_type = FORWARD_RECURRENCE_TIME;
      }
      else if (*pstr == "Mixture") {
        dist_type = NB_EVENT_MIXTURE;
      }
      else {
        genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                    "InterEvent or LengthBias or Backward or Forward or NbEvent or Mixture");
        return AMObj(AMObjType::ERROR);
      }
    }

    dist = ((Renewal*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , time);

    if (dist) {
      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV) ||
      (args[0].tag() == AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV) || (args[0].tag() == AMObjType::HIDDEN_SEMI_MARKOV) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV)) {
    string *pstr;
    bool status = true;
    process_distribution dist_type;
    DiscreteParametricModel *dist;
    StatError error;


    CHECKCONDVA((args.length() >= 3) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractDistribution"));

    if (args[1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                  args[1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[1].val.p;

      if (*pstr == "Observation") {
        dist_type = OBSERVATION;
      }
      else if (*pstr == "FirstOccurrence") {
        dist_type = FIRST_OCCURRENCE;
      }
      else if (*pstr == "Recurrence") {
        dist_type = RECURRENCE_TIME;
      }
      else if (*pstr == "Sojourn") {
        dist_type = SOJOURN_TIME;
      }
      else if (*pstr == "NbRun") {
        dist_type = NB_RUN;
      }
      else if (*pstr == "NbOccurrence") {
        dist_type = NB_OCCURRENCE;
      }
      else if (*pstr != "Forward") {
        status = false;
        genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                    "Observation or FirstOccurrence or Recurrence or Sojourn or Forward or NbRun or NbOccurrence");
      }
    }

    if (*pstr == "Forward") {
      int nb_required = 3;
      process_distribution histo_type = FINAL_RUN;


      CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractDistribution"));

      // arguments obligatoires

      if ((args[0].tag() != AMObjType::SEMI_MARKOV) && (args[0].tag() != AMObjType::HIDDEN_SEMI_MARKOV)) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 1 ,
                    args[0].tag.string().data() , "SEMI-MARKOV or HIDDEN_SEMI-MARKOV");
      }

      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 3 ,
                    args[2].tag.string().data() , "INT");
      }

      // argument optionnel

      if (args.length() == nb_required + 2) {
        if (args[nb_required].tag() != AMObjType::OPTION) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , nb_required + 1 ,
                      args[nb_required].tag.string().data() , "OPTION");
        }
        else {
          if (*((AMString*)args[nb_required].val.p) != "FrequencyDistributionType") {
            status = false;
            genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "ExtractDistribution" , nb_required + 1 ,
                        "FrequencyDistributionType");
          }
        }

        if (args[nb_required + 1].tag() != AMObjType::STRING) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , nb_required + 1 ,
                      args[nb_required + 1].tag.string().data() , "STRING");
        }
        else {
          pstr = (AMString*)args[nb_required + 1].val.p;
          if (*pstr == "InitialRun") {
            histo_type = INITIAL_RUN;
          }
          else if (*pstr == "FinalRun") {
            histo_type = FINAL_RUN;
          }
          else {
            status = false;
            genAMLError(ERRORMSG(ARG_sds) , "ExtractDistribution" , nb_required + 1 ,
                        "InitialRun or FinalRun");
          }
        }
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      switch (args[0].tag()) {
      case AMObjType::SEMI_MARKOV :
        dist = ((SemiMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[2].val.i , histo_type);
        break;
      case AMObjType::HIDDEN_SEMI_MARKOV :
        dist = ((HiddenSemiMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[2].val.i , histo_type);
        break;
      }
    }

    else {
      int variable , value;


      CHECKCONDVA((args.length() == 3) || (args.length() == 4) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractDistribution"));

      if ((dist_type == OBSERVATION) && (args.length() != 4)) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 4);
      }

      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 3 ,
                    args[2].tag.string().data() , "INT");
      }
      else if ((args.length() == 4) && (args[2].val.i <= 0)) {
        status = false;
        genAMLError(ERRORMSG(VARIABLE_INDEX_1_sd) , "ExtractDistribution" , 3);
      }

      if ((args.length() == 4) && (args[3].tag() != AMObjType::INTEGER)) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 4 ,
                    args[3].tag.string().data() , "INT");
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      switch (args.length()) {
      case 3 :
        variable = 0;
        value = args[2].val.i;
        break;
      case 4 :
        variable = args[2].val.i;
        value = args[3].val.i;
        break;
      }

      switch (args[0].tag()) {
      case AMObjType::VARIABLE_ORDER_MARKOV :
        dist = ((VariableOrderMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , variable , value);
        break;
      case AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV :
        dist = ((HiddenVariableOrderMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , variable , value);
        break;
      case AMObjType::SEMI_MARKOV :
        dist = ((SemiMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , variable , value);
        break;
      case AMObjType::HIDDEN_SEMI_MARKOV :
        dist = ((HiddenSemiMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , variable , value);
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV :
        dist = ((NonhomogeneousMarkov*)((STAT_model*)args[0].val.p)->pt)->extract(error , dist_type , value);
        break;
      }
    }

    if (dist) {
      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  genAMLError(ERRORMSG(STAT_MODEL_TYPE_DISTRIBUTION_sds) , "ExtractDistribution" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Extraction d'une loi empirique elementaire.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractFrequencyDistribution(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

  if (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) {
    string *pstr;
    DiscreteDistributionData *histo;
    const DiscreteMixtureData *mixture_data;
    StatError error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

    mixture_data = (DiscreteMixtureData*)((STAT_model*)args[0].val.p)->pt;

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Weight") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

      histo = new DiscreteDistributionData(*(mixture_data->get_weight()) ,
                                           mixture_data->get_mixture()->get_weight());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    if (*pstr == "Component") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 3 ,
                              args[2].tag.string().data() , "INT"));

      histo = mixture_data->extract(error , args[2].val.i);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Mixture") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

      histo = new DiscreteDistributionData(*((FrequencyDistribution*)mixture_data) ,
                                           (Distribution*)mixture_data->get_mixture());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                "Weight or Component or Mixture");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::CONVOLUTION_DATA) {
    string *pstr;
    DiscreteDistributionData *histo;
    const ConvolutionData *convolution_data;
    StatError error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

    convolution_data = (ConvolutionData*)((STAT_model*)args[0].val.p)->pt;

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Elementary") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 3 ,
                            args[2].tag.string().data() , "INT"));

      histo = convolution_data->extract(error , args[2].val.i);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Convolution") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

      histo = new DiscreteDistributionData(*((FrequencyDistribution*)convolution_data) ,
                                           (Distribution*)convolution_data->get_convolution());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                "Elementary or Convolution");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::COMPOUND_DATA) {
    string *pstr;
    compound_distribution type;
    DiscreteDistributionData *histo;
    StatError error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Sum") {
      type = SUM;
    }
    else if (*pstr == "Elementary") {
      type = ELEMENTARY;
    }
    else if (*pstr == "Compound") {
      type = COMPOUND;
    }
    else {
      genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                  "Sum or Elementary or Compound");
      return AMObj(AMObjType::ERROR);
    }

    histo = ((CompoundData*)((STAT_model*)args[0].val.p)->pt)->extract(error , type);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    int nb_variable , variable;
    DiscreteDistributionData *histo;
    const Vectors *vec;
    StatError error;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      vec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      vec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    nb_variable = vec->get_nb_variable();

    if (nb_variable == 1) {
      CHECKCONDVA(args.length() == 1 ,
                  genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ExtractFrequencyDistribution"));
      variable = 1;
    }

    else {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

      if (args[1].tag() != AMObjType::INTEGER) {
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                    args[1].tag.string().data() , "INT");
        return AMObj(AMObjType::ERROR);
      }
      else {
        variable = args[1].val.i;
      }
    }

    histo = vec->extract(error , variable);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::TIME_EVENTS) || (args[0].tag() == AMObjType::RENEWAL_DATA)) {
    string *pstr;
    renewal_distribution histo_type;
    DiscreteDistributionData *histo;
    const TimeEvents *time_events;
    const RenewalData *renewal_data;
    StatError error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    time_events = (TimeEvents*)((STAT_model*)args[0].val.p)->pt;

    if (args[0].tag() == AMObjType::RENEWAL_DATA) {
      renewal_data = (RenewalData*)((STAT_model*)args[0].val.p)->pt;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "NbEvent") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 3 ,
                              args[2].tag.string().data() , "INT"));

      switch (args[0].tag()) {
      case AMObjType::TIME_EVENTS :
        histo = time_events->extract(error , NB_EVENT , args[2].val.i);
        break;
      case AMObjType::RENEWAL_DATA :
        histo = renewal_data->extract(error , NB_EVENT , args[2].val.i);
        break;
      }

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

    if (*pstr == "Mixture") {
      switch (args[0].tag()) {
      case AMObjType::TIME_EVENTS :
        histo = time_events->extract(error , NB_EVENT_MIXTURE);
        break;
      case AMObjType::RENEWAL_DATA :
        histo = renewal_data->extract(error , NB_EVENT_MIXTURE);
        break;
      }

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    if (*pstr == "ObservationTime") {
      histo = new DiscreteDistributionData(*(time_events->get_htime()));

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    if (args[0].tag() == AMObjType::TIME_EVENTS) {
      genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                  "ObservationTime or NbEvent or Mixture");
      return AMObj(AMObjType::ERROR);
    }

    CHECKCONDVA(args[0].tag() == AMObjType::RENEWAL_DATA ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 1 ,
                            args[0].tag.string().data() , "RENEWAL_DATA"));

    if (*pstr == "InterEvent") {
      histo_type = INTER_EVENT;
    }
    else if (*pstr == "Within") {
      histo_type = WITHIN_OBSERVATION_PERIOD;
    }
    else if (*pstr == "LengthBias") {
      histo_type = LENGTH_BIAS;
    }
    else if (*pstr == "Backward") {
      histo_type = BACKWARD_RECURRENCE_TIME;
    }
    else if (*pstr == "Forward") {
      histo_type = FORWARD_RECURRENCE_TIME;
    }
    else {
      genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                  "ObservationTime or InterEvent or Within or LengthBias or Backward or Forward or NbEvent or Mixture");
      return AMObj(AMObjType::ERROR);
    }

    histo = renewal_data->extract(error , histo_type);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    string *pstr;
    bool status = true;
    int nb_variable , variable , offset , value;
    process_distribution histo_type;
    DiscreteDistributionData *histo;
    const Sequences *seq;
    const MarkovianSequences *markovian_seq;
    StatError error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      seq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Length") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

      histo = new DiscreteDistributionData(*(seq->get_length_distribution()));

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }

    if (*pstr == "IndexParameter") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));

//      if (seq->get_index_parameter_type() == IMPLICIT_TYPE) {
      if (!seq->get_index_parameter_distribution()) {
        genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NOT_BUILT_s) , "ExtractFrequencyDistribution");
        return AMObj(AMObjType::ERROR);
      }

      else {
        histo = new DiscreteDistributionData(*(seq->get_index_parameter_distribution()));
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
      }
    }

    if (*pstr == "Value") {
      nb_variable = seq->get_nb_variable();

      if (nb_variable == 1) {
        CHECKCONDVA(args.length() == 2 ,
                    genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 2));
        variable = 1;
      }

      else {
        CHECKCONDVA(args.length() == 3 ,
                    genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , 3));

        if (args[2].tag() != AMObjType::INTEGER) {
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 3 ,
                      args[2].tag.string().data() , "INT");
          return AMObj(AMObjType::ERROR);
        }
        else {
          variable = args[2].val.i;
        }
      }

      histo = seq->extract(error , variable);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
        return AMObj(AMObjType::ERROR);
      }
    }

    CHECKCONDVA(args.length() >= 3 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractFrequencyDistribution"));

    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      markovian_seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      markovian_seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      markovian_seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      markovian_seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    default :
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 1 , args[0].tag.string().data() ,
                  "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
      return AMObj(AMObjType::ERROR);
    }

    if (*pstr == "Observation") {
      histo_type = OBSERVATION;
    }
    else if (*pstr == "FirstOccurrence") {
      histo_type = FIRST_OCCURRENCE;
    }
    else if (*pstr == "Recurrence") {
      histo_type = RECURRENCE_TIME;
    }
    else if (*pstr == "Sojourn") {
      histo_type = SOJOURN_TIME;
    }
    else if (*pstr == "InitialRun") {
      histo_type = INITIAL_RUN;
    }
    else if (*pstr == "FinalRun") {
      histo_type = FINAL_RUN;
    }
    else if (*pstr == "NbRun") {
      histo_type = NB_RUN;
    }
    else if (*pstr == "NbOccurrence") {
      histo_type = NB_OCCURRENCE;
    }
    else {
      status = false;
      genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractFrequencyDistribution" , 2 ,
                  "Length or Observation or FirstOccurrence or Recurrence or Sojourn or InitialRun or FinalRun or NbRun or NbOccurrence");
    }

    if ((histo_type == OBSERVATION) && (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 1 ,
                  args[0].tag.string().data() , "VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    }

    nb_variable = markovian_seq->get_nb_variable();

    if (nb_variable == 1) {
      offset = 2;
      variable = 1;
    }

    else {
      offset = 3;

      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , 3 ,
                    args[2].tag.string().data() , "INT");
      }
      else {
        variable = args[2].val.i;
      }
    }

    CHECKCONDVA(args.length() == offset + 1 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractFrequencyDistribution" , offset + 1));

    if (args[offset].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractFrequencyDistribution" , offset + 1 ,
                  args[offset].tag.string().data() , "INT");
    }
    else {
      value = args[offset].val.i;
    }

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      histo = ((MarkovianSequences*)((STAT_model*)args[0].val.p)->pt)->extract(error , histo_type , variable , value);
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      histo = ((VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt)->extract(error , histo_type , variable , value);
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      histo = ((SemiMarkovData*)((STAT_model*)args[0].val.p)->pt)->extract(error , histo_type , variable , value);
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      histo = ((NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt)->extract(error , histo_type , value);
      break;
    }

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractFrequencyDistribution");
      return AMObj(AMObjType::ERROR);
    }
  }

  genAMLError(ERRORMSG(STAT_DATA_TYPE_FREQUENCY_DISTRIBUTION_sds) , "ExtractFrequencyDistribution" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Pour compatibilite ascendante : Extraction d'une loi empirique elementaire.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractHistogram(const AMObjVector &args)

{
  return STAT_ExtractFrequencyDistribution(args);
}


/*--------------------------------------------------------------*
 *
 *  Extraction de mesures globales (longueur, temps avant la 1ere occurrence
 *  d'une valeur, nombre de series/d'occurrences d'une valeur) par sequence.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractVectors(const AMObjVector &args)

{
  string *pstr;
  bool status = true;
  int nb_variable , variable , offset , value;
  sequence_pattern pattern; 
  Vectors *vec;
  const Sequences *seq;
  StatError error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractVectors"));

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    seq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractVectors" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractVectors" , 2 ,
                          args[1].tag.string().data() , "STRING"));

  pstr = (AMString*)args[1].val.p;

  if (*pstr == "Length") {
    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractVectors" , 2));

    vec = seq->extract_vectors(error , LENGTH_PATTERN);

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractVectors");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (*pstr == "Cumul") {
    pattern = SEQUENCE_CUMUL;
  }
  else if (*pstr == "Mean") {
    pattern = SEQUENCE_MEAN;
  }
  else if (*pstr == "FirstOccurrence") {
    pattern = FIRST_OCCURRENCE_PATTERN;
  }
  else if (*pstr == "Sojourn") {
    pattern = SOJOURN_TIME_PATTERN;
  }
  else if (*pstr == "NbRun") {
    pattern = NB_RUN_PATTERN;
  }
  else if (*pstr == "NbOccurrence") {
    pattern = NB_OCCURRENCE_PATTERN;
  }
  else {
    status = false;
    genAMLError(ERRORMSG(FREQUENCY_DISTRIBUTION_NAME_sds) , "ExtractVectors" , 2 ,
                "Length or Cumul or Mean or FirstOccurrence or NbRun or NbOccurrence");
  }

  nb_variable = seq->get_nb_variable();

  if (nb_variable == 1) {
    offset = 2;
    variable = 1;
  }

  else {
    offset = 3;

    if (args[2].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractVectors" , 3 ,
                    args[2].tag.string().data() , "INT");
    }
    else {
      variable = args[2].val.i;
    }
  }

  if ((*pstr == "Cumul") || (*pstr == "Mean")) {
    CHECKCONDVA(args.length() == offset ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractVectors" , offset));

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    vec = seq->extract_vectors(error , pattern , variable);

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractVectors");
      return AMObj(AMObjType::ERROR);
    }
  }

  CHECKCONDVA(args.length() == offset + 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractVectors" , offset + 1));

  if (args[offset].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractVectors" , offset + 1 ,
                args[offset].tag.string().data() , "INT");
  }
  else {
    value = args[offset].val.i;
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  vec = seq->extract_vectors(error , pattern , variable , value);

  if (vec) {
    STAT_model* model = new STAT_model(vec);
    return AMObj(AMObjType::VECTORS , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractVectors");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction des parties "donnees" d'un objet de type "modele".
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractData(const AMObjVector &args)

{
  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ExtractData"));

  switch (args[0].tag()) {

  case AMObjType::DISCRETE_MIXTURE : {
    DiscreteMixtureData *mixt_histo;
    StatError error;


    mixt_histo = ((DiscreteMixture*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (mixt_histo) {
      STAT_model* model = new STAT_model(mixt_histo);
      return AMObj(AMObjType::DISCRETE_MIXTURE_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::CONVOLUTION : {
    ConvolutionData *convol_histo;
    StatError error;


    convol_histo = ((Convolution*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (convol_histo) {
      STAT_model* model = new STAT_model(convol_histo);
      return AMObj(AMObjType::CONVOLUTION_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::COMPOUND : {
    CompoundData *compound_histo;
    StatError error;


    compound_histo = ((Compound*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (compound_histo) {
      STAT_model* model = new STAT_model(compound_histo);
      return AMObj(AMObjType::COMPOUND_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::MIXTURE : {
    MixtureData *vec;
    StatError error;


    vec = ((Mixture*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::MIXTURE_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV : {
    VariableOrderMarkovData *seq;
    StatError error;


    seq = ((HiddenVariableOrderMarkov*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::VARIABLE_ORDER_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::HIDDEN_SEMI_MARKOV : {
    SemiMarkovData *seq;
    StatError error;


    seq = ((HiddenSemiMarkov*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEMI_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  default : {
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ExtractData" , args[0].tag.string().data() ,
                "DISCRETE_MIXTURE or CONVOLUTION or COMPOUND or MIXTURE or HIDDEN_VARIABLE_ORDER_MARKOV or HIDDEN_SEMI-MARKOV");
    return AMObj(AMObjType::ERROR);
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Fusion de structures de donnees de meme type ou de correlogrammes.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Merge(const AMObjVector &args)

{
  bool status = true;
  int i;
  int nb_sample = args.length();


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Merge"));

  if ((args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION) || (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    const FrequencyDistribution **phisto;
    vector<FrequencyDistribution> ihisto;
    DiscreteDistributionData *histo;


    phisto = new const FrequencyDistribution*[nb_sample];
    for (i = 1;i < nb_sample;i++) {
      phisto[i] = NULL;
    }

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::FREQUENCY_DISTRIBUTION :
        phisto[i] = (FrequencyDistribution*)((DiscreteDistributionData*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::DISCRETE_MIXTURE_DATA :
        phisto[i] = (FrequencyDistribution*)((DiscreteMixtureData*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::CONVOLUTION_DATA :
        phisto[i] = (FrequencyDistribution*)((ConvolutionData*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::COMPOUND_DATA :
        phisto[i] = (FrequencyDistribution*)((CompoundData*)((STAT_model*)args[i].val.p)->pt);
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 , args[i].tag.string().data() ,
                    "FREQUENCY_DISTRIBUTION or DISCRETE_MIXTURE_DATA or CONVOLUTION_DATA or COMPOUND_DATA");
        break;
      }

      if ((i > 0) && (phisto[i])) {
        ihisto.push_back(*phisto[i]);
      }
    }

    if (status) {
//      histo = new DiscreteDistributionData(nb_sample , phisto);
      histo = phisto[0]->merge(nb_sample - 1 , ihisto);
    }

    delete [] phisto;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(histo);
    return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors **pvec;
    vector<Vectors> ivec;
    Vectors *vec;
    StatError error;


    pvec = new const Vectors*[nb_sample];
    for (i = 1;i < nb_sample;i++) {
      pvec[i] = NULL;
    }

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::VECTORS :
        pvec[i] = (Vectors*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MIXTURE_DATA :
        pvec[i] = (MixtureData*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "VECTORS or MIXTURE_DATA");
      }

      if ((i > 0) && (pvec[i])) {
        ivec.push_back(*pvec[i]);
      }
    }

    if (status) {
//      vec = pvec[0]->merge(error , nb_sample - 1 , pvec + 1);
      vec = pvec[0]->merge(error , nb_sample - 1 , ivec);
    }

    delete [] pvec;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Merge");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::TIME_EVENTS) {
    const TimeEvents **ptimev;
    vector<TimeEvents> itimev;
    TimeEvents *timev;


    ptimev = new const TimeEvents*[nb_sample];
    for (i = 1;i < nb_sample;i++) {
      ptimev[i] = NULL;
    }

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::TIME_EVENTS :
        ptimev[i] = (TimeEvents*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::RENEWAL_DATA :
        ptimev[i] = (RenewalData*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "TIME_EVENTS or RENEWAL_DATA");
        break;
      }

      if ((i > 0) && (ptimev[i])) {
        itimev.push_back(*ptimev[i]);
      }
    }

    if (status) {
//      timev = new TimeEvents(nb_sample , ptimev);
      timev = ptimev[0]->merge(nb_sample - 1 , itimev);
    }

    delete [] ptimev;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(timev);
    return AMObj(AMObjType::TIME_EVENTS , model);
  }

  if (args[0].tag() == AMObjType::RENEWAL_DATA) {
    const RenewalData **ptimev;
    RenewalData *timev;
    vector<RenewalData> itimev;
    StatError error;


    ptimev = new const RenewalData*[nb_sample];
    for (i = 1;i < nb_sample;i++) {
      ptimev[i] = NULL;
    }

    for (i = 0;i < nb_sample;i++) {
      if (args[i].tag() != AMObjType::RENEWAL_DATA) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "RENEWAL_DATA");
      }
      else {
        ptimev[i] = (RenewalData*)((STAT_model*)args[i].val.p)->pt;
      }

      if ((i > 0) && (ptimev[i])) {
        itimev.push_back(*ptimev[i]);
      }
    }

    if (status) {
//      timev = ptimev[0]->merge(error , nb_sample - 1 , ptimev + 1);
      timev = ptimev[0]->merge(error , nb_sample - 1 , itimev);
    }

    delete [] ptimev;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (timev) {
      STAT_model* model = new STAT_model(timev);
      return AMObj(AMObjType::RENEWAL_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Merge");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences **pseq;
    vector<Sequences> iseq;
    Sequences *seq;
    const MarkovianSequences **pmarkovian_seq;
    vector<MarkovianSequences> imarkovian_seq;
    MarkovianSequences *markovian_seq;
    StatError error;


    pseq = new const Sequences*[nb_sample];
    for (i = 1;i < nb_sample;i++) {
      pseq[i] = NULL;
    }

    pmarkovian_seq = new const MarkovianSequences*[nb_sample];
    for (i = 0;i < nb_sample;i++) {
      pmarkovian_seq[i] = NULL;
    }

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::SEQUENCES :
        pseq[i] = (Sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MARKOVIAN_SEQUENCES :
        pseq[i] = (MarkovianSequences*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (MarkovianSequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
        pseq[i] = (VariableOrderMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (VariableOrderMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::SEMI_MARKOV_DATA :
        pseq[i] = (SemiMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (SemiMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
        pseq[i] = (NonhomogeneousMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (NonhomogeneousMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 , args[i].tag.string().data() ,
                    "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
        break;
      }

      if (i > 0) {
        if (pseq[i]) {
          iseq.push_back(*pseq[i]);
        }
        if (pmarkovian_seq[i]) {
          imarkovian_seq.push_back(*pmarkovian_seq[i]);
        }
      }
    }

    if (!status) {
      delete [] pseq;
      delete [] pmarkovian_seq;
      return AMObj(AMObjType::ERROR);
    }

    for (i = 0;i < nb_sample;i++) {
      if (!pmarkovian_seq[i]) {
        break;
      }
    }

    if (i < nb_sample) {
//      seq = pseq[0]->merge(error , nb_sample - 1 , pseq + 1);
      seq = pseq[0]->merge(error , nb_sample - 1 , iseq);
      delete [] pseq;
      delete [] pmarkovian_seq;

      if (seq) {
        markovian_seq = seq->markovian_sequences(error);
        if (markovian_seq) {
          delete seq;
          STAT_model* model = new STAT_model(markovian_seq);
          return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
        }
        else {
          AMLOUTPUT << "\n";
          error.ascii_write(AMLOUTPUT , WARNING);
          STAT_model* model = new STAT_model(seq);
          return AMObj(AMObjType::SEQUENCES , model);
        }
      }

      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Merge");
        return AMObj(AMObjType::ERROR);
      }
    }

    else {
//      markovian_seq = pmarkovian_seq[0]->merge(error , nb_sample - 1 , pmarkovian_seq + 1);
      markovian_seq = pmarkovian_seq[0]->merge(error , nb_sample - 1 , imarkovian_seq);
      delete [] pseq;
      delete [] pmarkovian_seq;

      if (markovian_seq) {
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Merge");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

  if (args[0].tag() == AMObjType::CORRELATION) {
    const Correlation **pcorrel;
    Correlation *correl;
    StatError error;


    pcorrel = new const Correlation*[nb_sample];
    pcorrel[0] = (Correlation*)((STAT_model*)args[0].val.p)->pt;

    for (i = 1;i < nb_sample;i++) {
      if (args[i].tag() != AMObjType::CORRELATION) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "CORRELATION");
      }
      else {
        pcorrel[i] = (Correlation*)((STAT_model*)args[i].val.p)->pt;
      }
    }

    if (status) {
      correl = pcorrel[0]->merge(error , nb_sample - 1 , pcorrel + 1);
    }

    delete [] pcorrel;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (correl) {
      STAT_model* model = new STAT_model(correl);
      return AMObj(AMObjType::CORRELATION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Merge");
      return AMObj(AMObjType::ERROR);
    }
  }

  genAMLError(ERRORMSG(STAT_DATA_TYPE_CORRELATION_sds) , "Merge" , 1 , args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Translation des valeurs.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Shift(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Shift"));

  if ((args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION) || (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    bool status = true;
    const FrequencyDistribution *ihisto;
    DiscreteDistributionData *histo;
    StatError error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Shift" , 2));

    switch (args[0].tag()) {
    case AMObjType::FREQUENCY_DISTRIBUTION :
      ihisto = (FrequencyDistribution*)((DiscreteDistributionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::DISCRETE_MIXTURE_DATA :
      ihisto = (FrequencyDistribution*)((DiscreteMixtureData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (FrequencyDistribution*)((ConvolutionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (FrequencyDistribution*)((CompoundData*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Shift" , 2 ,
                  args[1].tag.string().data() , "INT");
    }

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    histo = ihisto->shift(error , args[1].val.i);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Shift");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    bool status = true;
    int nb_variable , variable , offset;
    const Vectors *ivec;
    const Sequences *iseq;
    StatError error;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    }

    if (nb_variable == 1) {
      offset = 1;
      variable = 1;
    }

    else {
      offset = 2;

      if (args[1].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Shift" , 2 ,
                    args[1].tag.string().data() , "INT");
      }
      else {
        variable = args[1].val.i;
      }
    }

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    CHECKCONDVA(args.length() == offset + 1 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Shift" , offset + 1));

    if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
      Vectors *vec;


      if (args[offset].tag() == AMObjType::INTEGER) {
        vec = ivec->shift(error , variable , args[offset].val.i);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        vec = ivec->shift(error , variable , args[offset].val.r);
      }
      else {
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Shift" , offset + 1 ,
                    args[offset].tag.string().data() , "INT or REAL");
        return AMObj(AMObjType::ERROR);
      }

      if (vec) {
        STAT_model* model = new STAT_model(vec);
        return AMObj(AMObjType::VECTORS , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Shift");
        return AMObj(AMObjType::ERROR);
      }
    }

    if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
        (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      Sequences *seq;
      MarkovianSequences *markovian_seq;


      if (args[offset].tag() == AMObjType::INTEGER) {
        seq = iseq->shift(error , variable , args[offset].val.i);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        seq = iseq->shift(error , variable , args[offset].val.r);
      }
      else {
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Shift" , offset + 1 ,
                    args[offset].tag.string().data() , "INT or REAL");
        return AMObj(AMObjType::ERROR);
      }

      if (seq) {
        markovian_seq = seq->markovian_sequences(error);
        if (markovian_seq) {
          delete seq;
          STAT_model* model = new STAT_model(markovian_seq);
          return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
        }
        else {
          AMLOUTPUT << "\n";
          error.ascii_write(AMLOUTPUT , WARNING);
          STAT_model* model = new STAT_model(seq);
          return AMObj(AMObjType::SEQUENCES , model);
        }
      }

      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Shift");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

//  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Shift" , 1 ,
  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_sds) , "Shift" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Seuillage des valeurs d'une variable.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ThresholdingData(const AMObjVector &args)

{
  string *pstr;
  bool status = true;
  int nb_variable , variable , offset;
  threshold_direction mode;
  const Vectors *ivec;
  const Sequences *iseq;
  StatError error;


  CHECKCONDVA(args.length() >= 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Thresholding"));

  switch (args[0].tag()) {
  case AMObjType::VECTORS :
    ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = ivec->get_nb_variable();
    break;
  case AMObjType::MIXTURE_DATA :
    ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = ivec->get_nb_variable();
    break;
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Thresholding" , 1 , args[0].tag.string().data() ,
                "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Thresholding" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  CHECKCONDVA(args.length() == offset + 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Thresholding" , offset + 2));

  if (args[offset + 1].tag() != AMObjType::STRING) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Thresholding" , offset + 2 ,
                args[offset + 1].tag.string().data() , "STRING");
  }
  else {
    pstr = (AMString*)args[offset + 1].val.p;
    if (*pstr == "Above") {
      mode = ABOVE;
    }
    else if (*pstr == "Below") {
      mode = BELOW;
    }
    else {
      status = false;
      genAMLError(ERRORMSG(ARG_sds) , "Thresholding" , offset + 2 , "Above or Below");
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    Vectors *vec;


    if (args[offset].tag() == AMObjType::INTEGER) {
      vec = ivec->thresholding(error , variable , args[offset].val.i , mode);
    }
    else if (args[offset].tag() == AMObjType::REAL) {
      vec = ivec->thresholding(error , variable , args[offset].val.r , mode);
    }
    else {
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Thresholding" , offset + 1 ,
                  args[offset].tag.string().data() , "INT or REAL");
      return AMObj(AMObjType::ERROR);
    }

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Thresholding");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    if (args[offset].tag() == AMObjType::INTEGER) {
      seq = iseq->thresholding(error , variable , args[offset].val.i , mode);
    }
    else if (args[offset].tag() == AMObjType::REAL) {
      seq = iseq->thresholding(error , variable , args[offset].val.r , mode);
    }
    else {
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Thresholding" , offset + 1 ,
                  args[offset].tag.string().data() , "INT or REAL");
      return AMObj(AMObjType::ERROR);
    }

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Thresholding");
      return AMObj(AMObjType::ERROR);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Groupement des classes.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Cluster(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cluster"));

  CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 2 ,
                          args[1].tag.string().data() , "STRING"));

  if ((args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION) || (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    string *pstr;
    bool status = true;
    int nb_required = 3;
    const FrequencyDistribution *ihisto;
    DiscreteDistributionData *histo;
    StatError error;


    switch (args[0].tag()) {
    case AMObjType::FREQUENCY_DISTRIBUTION :
      ihisto = (FrequencyDistribution*)((DiscreteDistributionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::DISCRETE_MIXTURE_DATA :
      ihisto = (FrequencyDistribution*)((DiscreteMixtureData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (FrequencyDistribution*)((ConvolutionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (FrequencyDistribution*)((CompoundData*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Step") {
      rounding mode = FLOOR;


      CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cluster"));

      // argument obligatoire

      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 3 ,
                    args[2].tag.string().data() , "INT");
      }

      // argument optionnel

      if (args.length() == nb_required + 2) {
        if (args[nb_required].tag() != AMObjType::OPTION) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                      args[nb_required].tag.string().data() , "OPTION");
        }
        else {
          if (*((AMString*)args[nb_required].val.p) != "Mode") {
            status = false;
            genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cluster" , nb_required + 1 , "Mode");
          }
        }

        if (args[nb_required + 1].tag() != AMObjType::STRING) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                      args[nb_required + 1].tag.string().data() , "STRING");
        }
        else {
          pstr = (AMString*)args[nb_required + 1].val.p;
          if (*pstr == "Floor") {
            mode = FLOOR;
          }
          else if (*pstr == "Round") {
            mode = ROUND;
          }
          else if (*pstr == "Ceil") {
            mode = CEIL;
          }
          else {
            status = false;
            genAMLError(ERRORMSG(ROUND_MODE_sds) , "Cluster" , nb_required + 1 ,
                        "Floor or Round or Ceil");
          }
        }
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      histo = ihisto->cluster(error , args[2].val.i , mode);
    }

    else if (*pstr == "Information") {
      double ratio;


      CHECKCONDVA(args.length() == nb_required ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , 3));

      switch (args[2].tag()) {
      case AMObjType::INTEGER :
        ratio = args[2].val.i;
        break;
      case AMObjType::REAL :
        ratio = args[2].val.r;
        break;
      default :
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 3 ,
                    args[2].tag.string().data() , "INT or REAL");
        return AMObj(AMObjType::ERROR);
      }

      histo = ihisto->cluster(error , ratio , &AMLOUTPUT);
    }

    else if (*pstr == "Limit") {
      int i;
      int nb_class = I_DEFAULT , *limit = NULL;
      vector<int> vec_limit;


      CHECKCONDVA(args.length() == nb_required ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , 3));

      if (args[2].tag() != AMObjType::ARRAY) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 3 ,
                    args[2].tag.string().data() , "ARRAY");
      }
      else {
        limit = buildIntArray(args , 2 , "Cluster" , 3 , nb_class);
        if (!limit) {
          status = false;
        }
        else {
          for (i = 0;i < nb_class;i++) {
            vec_limit.push_back(limit[i]);
          }
        }
      }

      if (status) {
//        histo = ihisto->cluster(error , nb_class + 1 , limit);
        histo = ihisto->cluster(error , nb_class + 1 , vec_limit);
      }

      delete [] limit;

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }
    }

    else {
      genAMLError(ERRORMSG(CLUSTERING_METHOD_sds) , "Cluster" , 2 ,
                  "Step or Information or Limit");
      return AMObj(AMObjType::ERROR);
    }

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Cluster");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    string *pstr;
    bool status = true;
    int i;
    int nb_required , nb_variable , variable , offset , nb_class = I_DEFAULT , *int_limit = NULL;
    vector<int> vec_int_limit;
    rounding  mode = FLOOR;
    double *real_limit = NULL;
    vector<double> vec_real_limit;
    const Vectors *ivec;
    const Sequences *iseq;
    const MarkovianSequences *imarkovian_seq;
    StatError error;


    // arguments obligatoires

    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      imarkovian_seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      imarkovian_seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      imarkovian_seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      imarkovian_seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    }

    if (nb_variable == 1) {
      offset = 2;
      variable = 1;
    }

    else {
      offset = 3;

      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 3 ,
                    args[2].tag.string().data() , "INT");
      }
      else {
        variable = args[2].val.i;
      }
    }

    nb_required = offset + 1;

    CHECKCONDVA(args.length() >= nb_required ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cluster"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Step") {
      if (args[offset].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , offset + 1 ,
                    args[offset].tag.string().data() , "INT");
      }
    }

    else if (*pstr == "Limit") {
      if (args[offset].tag() == AMObjType::ARRAY) {
        Array* parray = (Array*)args[offset].val.p;

        if (parray->surfaceType() == AMObjType::INTEGER) {
          int_limit = buildIntArray(args , offset , "Cluster" , offset + 1 , nb_class);
          if (!int_limit) {
            status = false;
          }
          else {
            for (i = 0;i < nb_class;i++) {
              vec_int_limit.push_back(int_limit[i]);
            }
          }
        }

        else if (parray->surfaceType() == AMObjType::REAL) {
          real_limit = buildRealArray(args , offset , "Cluster" , offset + 1 , nb_class , false);
          if (!real_limit) {
            status = false;
          }
          else {
            for (i = 0;i < nb_class;i++) {
              vec_real_limit.push_back(real_limit[i]);
            }
          }
        }

        else {
          status = false;
          genAMLError(ERRORMSG(ARRAY_ELEMENT_TYPE_1_sdss) , "Cluster" , offset + 1 ,
                      (parray->surfaceType()).string().data() , "INT or REAL");
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , offset + 1 ,
                    args[offset].tag.string().data() , "ARRAY");
      }
    }

    else {
      status = false;
      genAMLError(ERRORMSG(CLUSTERING_METHOD_sds) , "Cluster" , 2 ,
                  "Step or Limit");
    }

    if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
      Vectors *vec;


      if (*pstr == "Step") {
        if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }

        // argument optionnel

        if (args.length() == nb_required + 2) {
          if (args[nb_required].tag() != AMObjType::OPTION) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required].tag.string().data() , "OPTION");
          }
          else {
            if (*((AMString*)args[nb_required].val.p) != "Mode") {
              status = false;
              genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cluster" , nb_required + 1 , "Mode");
            }
          }

          if (args[nb_required + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + 1].val.p;
            if (*pstr == "Floor") {
              mode = FLOOR;
            }
            else if (*pstr == "Round") {
              mode = ROUND;
            }
            else if (*pstr == "Ceil") {
              mode = CEIL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(ROUND_MODE_sds) , "Cluster" , nb_required + 1 ,
                          "Floor or Round or Ceil");
            }
          }
        }

        if (status) {
          vec = ivec->cluster(error , variable , args[offset].val.i , mode);
        }
      }

      else if (*pstr == "Limit") {
        if (args.length() != nb_required) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }
        
        if (status) {
          if (int_limit) {
//            vec = ivec->cluster(error , variable , nb_class + 1 , int_limit);
            vec = ivec->cluster(error , variable , nb_class + 1 , vec_int_limit);
          }
          else if (real_limit) {
//            vec = ivec->cluster(error , variable , nb_class + 1 , real_limit);
            vec = ivec->cluster(error , variable , nb_class + 1 , vec_real_limit);
          }
        }

        delete [] int_limit;
        delete [] real_limit;
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (vec) {
        STAT_model* model = new STAT_model(vec);
        return AMObj(AMObjType::VECTORS , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Cluster");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (args[0].tag() == AMObjType::SEQUENCES) {
      Sequences *seq;
      MarkovianSequences *markovian_seq;


      if (*pstr == "Step") {
        if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }

        // argument optionnel

        if (args.length() == nb_required + 2) {
          if (args[nb_required].tag() != AMObjType::OPTION) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required].tag.string().data() , "OPTION");
          }
          else {
            if (*((AMString*)args[nb_required].val.p) != "Mode") {
              status = false;
              genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cluster" , nb_required + 1 , "Mode");
            }
          }

          if (args[nb_required + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + 1].val.p;
            if (*pstr == "Floor") {
              mode = FLOOR;
            }
            else if (*pstr == "Round") {
              mode = ROUND;
            }
            else if (*pstr == "Ceil") {
              mode = CEIL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(ROUND_MODE_sds) , "Cluster" , nb_required + 1 ,
                          "Floor or Round or Ceil");
            }
          }
        }

        if (status) {
          seq = iseq->cluster(error , variable , args[offset].val.i , mode);
        }
      }

      else if (*pstr == "Limit") {
        if (args.length() != nb_required) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }

        if (status) {
          if (int_limit) {
//            seq = iseq->cluster(error , variable , nb_class + 1 , int_limit);
            seq = iseq->cluster(error , variable , nb_class + 1 , vec_int_limit);
          }
          else if (real_limit) {
//            seq = iseq->cluster(error , variable , nb_class + 1 , real_limit);
            seq = iseq->cluster(error , variable , nb_class + 1 , vec_real_limit);
          }
        }

        delete [] int_limit;
        delete [] real_limit;
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (seq) {
        markovian_seq = seq->markovian_sequences(error);
        if (markovian_seq) {
          delete seq;
          STAT_model* model = new STAT_model(markovian_seq);
          return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
        }
        else {
          AMLOUTPUT << "\n";
          error.ascii_write(AMLOUTPUT , WARNING);
          STAT_model* model = new STAT_model(seq);
          return AMObj(AMObjType::SEQUENCES , model);
        }
      }

      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Cluster");
        return AMObj(AMObjType::ERROR);
      }
    }

    if ((args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
        (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      bool add_variable = false;
      MarkovianSequences *markovian_seq;


      pstr = (AMString*)args[1].val.p;

      if (*pstr == "Step") {
        if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }

        // argument optionnel

        if (args.length() == nb_required + 2) {
          if (args[nb_required].tag() != AMObjType::OPTION) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required].tag.string().data() , "OPTION");
          }
          else {
            if (*((AMString*)args[nb_required].val.p) != "Mode") {
              status = false;
              genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cluster" , nb_required + 1 , "Mode");
            }
          }

          if (args[nb_required + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + 1].val.p;
            if (*pstr == "Floor") {
              mode = FLOOR;
            }
            else if (*pstr == "Round") {
              mode = ROUND;
            }
            else if (*pstr == "Ceil") {
              mode = CEIL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(ROUND_MODE_sds) , "Cluster" , nb_required + 1 ,
                          "Floor or Round or Ceil");
            }
          }
        }

        if (status) {
          markovian_seq = imarkovian_seq->cluster(error , variable , args[offset].val.i , mode);
        }
      }

      else if (*pstr == "Limit") {
        if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
          delete [] int_limit;
          delete [] real_limit;

          genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cluster");
          return AMObj(AMObjType::ERROR);
        }

        // argument optionnel

        if (args.length() == nb_required + 2) {
          if (args[nb_required].tag() != AMObjType::OPTION) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required].tag.string().data() , "OPTION");
          }
          else {
            if (*((AMString*)args[nb_required].val.p) != "AddVariable") {
              status = false;
              genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cluster" , nb_required + 1 , "AddVariable");
            }
          }

          if (args[nb_required + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , nb_required + 1 ,
                        args[nb_required + 1].tag.string().data() , "BOOL");
          }
          else {
            add_variable = args[nb_required + 1].val.b;
          }

          if (real_limit) {
            status = false;
            genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "Cluster" , "AddVariable");
          }
        }

        if (status) {
          if (int_limit) {
//            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
//                                                    int_limit , add_variable);
            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
                                                    vec_int_limit , add_variable);
          }
          else if (real_limit) {
//            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
//                                                    real_limit);
            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
                                                    vec_real_limit);
          }
        }

        delete [] int_limit;
        delete [] real_limit;
      }

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (markovian_seq) {
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Cluster");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

//  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Cluster" , 1 ,
  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_sds) , "Cluster" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des categories.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Transcode(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Transcode"));

  if ((args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION) || (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    bool status = true;
    int i;
    int nb_category , *category = NULL;
    vector<int> vec_category;
    const FrequencyDistribution *ihisto;
    DiscreteDistributionData *histo;
    StatError error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , 2));

    switch (args[0].tag()) {
    case AMObjType::FREQUENCY_DISTRIBUTION :
      ihisto = (FrequencyDistribution*)((DiscreteDistributionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::DISCRETE_MIXTURE_DATA :
      ihisto = (FrequencyDistribution*)((DiscreteMixtureData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (FrequencyDistribution*)((ConvolutionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (FrequencyDistribution*)((CompoundData*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    nb_category = ihisto->nb_value - ihisto->offset;

    if (args[1].tag() != AMObjType::ARRAY) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , 2 ,
                  args[1].tag.string().data() , "ARRAY");
    }
    else {
      category = buildIntArray(args , 1 , "Transcode" , 2 , nb_category);
      if (!category) {
        status = false;
      }
      else {
        for (i = 0;i < nb_category;i++) {
          vec_category.push_back(category[i]);
        }
      }
    }

    if (status) {
//      histo = ihisto->transcode(error , category);
      histo = ihisto->transcode(error , vec_category);
    }

    delete [] category;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Transcode");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    bool status = true;
    int i;
    int nb_required , nb_variable , variable = I_DEFAULT , offset ,
        nb_category = I_DEFAULT , *category = NULL;
    vector<int> vec_category;
    const Vectors *ivec = NULL;
    const Sequences *iseq = NULL;
    const MarkovianSequences *imarkovian_seq = NULL;
    StatError error;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      imarkovian_seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      imarkovian_seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      imarkovian_seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      imarkovian_seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    }

    if (nb_variable == 1) {
      offset = 1;
      variable = 1;
    }

    else {
      offset = 2;

      if (args[1].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , 2 ,
                    args[1].tag.string().data() , "INT");
      }
      else {
         variable = args[1].val.i;
         CHECKCONDVA((variable >= 1) && (variable <= nb_variable) ,
                     genAMLError(ERRORMSG(VARIABLE_INDEX_2_sd) , "Transcode" , 2));
      }
    }

    nb_required = offset + 1;

    CHECKCONDVA(args.length() >= nb_required ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Transcode"));

    if (ivec) {
      nb_category = (int)(ivec->get_max_value(variable - 1) - ivec->get_min_value(variable - 1)) + 1;
    }
    else if (iseq) {
      nb_category = (int)(iseq->get_max_value(variable - 1) - iseq->get_min_value(variable - 1)) + 1;
    }
    else {
      nb_category = (int)(imarkovian_seq->get_max_value(variable - 1) - imarkovian_seq->get_min_value(variable - 1)) + 1;
    }

    if (args[offset].tag() != AMObjType::ARRAY) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , offset + 1 ,
                  args[offset].tag.string().data() , "ARRAY");
    }
    else {
      category = buildIntArray(args , offset , "Transcode" , offset + 1 , nb_category);
      if (!category) {
        status = false;
      }
      else {
        for (i = 0;i < nb_category;i++) {
          vec_category.push_back(category[i]);
        }
      }
    }

    if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
      Vectors *vec;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , nb_required);
      }

      if (status) {
//        vec = ivec->transcode(error , variable , category);
        vec = ivec->transcode(error , variable , vec_category);
      }

      delete [] category;

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (vec) {
        STAT_model* model = new STAT_model(vec);
        return AMObj(AMObjType::VECTORS , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Transcode");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (args[0].tag() == AMObjType::SEQUENCES) {
      Sequences *seq;
      MarkovianSequences *markovian_seq;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , nb_required);
      }

      if (status) {
//        seq = iseq->transcode(error , variable , category);
        seq = iseq->transcode(error , variable , vec_category);
      }

      delete [] category;

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (seq) {
        markovian_seq = seq->markovian_sequences(error);
        if (markovian_seq) {
          delete seq;
          STAT_model* model = new STAT_model(markovian_seq);
          return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
        }
        else {
          AMLOUTPUT << "\n";
          error.ascii_write(AMLOUTPUT , WARNING);
          STAT_model* model = new STAT_model(seq);
          return AMObj(AMObjType::SEQUENCES , model);
        }
      }

      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Transcode");
        return AMObj(AMObjType::ERROR);
      }
    }

    if ((args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
        (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      bool add_variable = false;
      MarkovianSequences *markovian_seq;


      if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
        delete [] category;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Transcode");
        return AMObj(AMObjType::ERROR);
      }

      // argument optionnel

      if (args.length() == nb_required + 2) {
        if (args[nb_required].tag() != AMObjType::OPTION) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , nb_required + 1 ,
                      args[nb_required].tag.string().data() , "OPTION");
        }
        else {
          if (*((AMString*)args[nb_required].val.p) != "AddVariable") {
            status = false;
            genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Transcode" , nb_required + 1 ,
                        "AddVariable");
          }
        }

        if (args[nb_required + 1].tag() != AMObjType::BOOL) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , nb_required + 1 ,
                      args[nb_required + 1].tag.string().data() , "BOOL");
        }
        else {
          add_variable = args[nb_required + 1].val.b;
        }
      }

      if (status) {
//        markovian_seq = imarkovian_seq->transcode(error , variable , category , add_variable);
        markovian_seq = imarkovian_seq->transcode(error , variable , vec_category , add_variable);
      }

      delete [] category;

      if (!status) {
        return AMObj(AMObjType::ERROR);
      }

      if (markovian_seq) {
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "Transcode");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

//  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Transcode" , 1 ,
  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_sds) , "Transcode" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Selection sur les valeurs prises par une variable.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ValueSelect(const AMObjVector &args)

{
  string *pstr;
  bool status = true , keep = true;
  int nb_required , int_min_value , int_max_value;
  StatError error;


  nb_required = nb_required_computation(args);

  CHECKCONDVA(nb_required >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "ValueSelect" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[nb_required + 1].val.p;
      if (*pstr == "Keep") {
        keep = true;
      }
      else if (*pstr == "Reject") {
        keep = false;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(ARG_sds) , "ValueSelect" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if ((args[0].tag() == AMObjType::FREQUENCY_DISTRIBUTION) || (args[0].tag() == AMObjType::DISCRETE_MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    const FrequencyDistribution *ihisto;
    DiscreteDistributionData *histo;


    CHECKCONDVA((nb_required == 2) || (nb_required == 3) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

    switch (args[0].tag()) {
    case AMObjType::FREQUENCY_DISTRIBUTION :
      ihisto = (FrequencyDistribution*)((DiscreteDistributionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::DISCRETE_MIXTURE_DATA :
      ihisto = (FrequencyDistribution*)((DiscreteMixtureData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (FrequencyDistribution*)((ConvolutionData*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (FrequencyDistribution*)((CompoundData*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      int_min_value = args[1].val.i;
    }

    switch (nb_required) {

    case 2 : {
      int_max_value = int_min_value;
      break;
    }

    case 3 : {
      if (args[2].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , 3 ,
                    args[2].tag.string().data() , "INT");
      }
      else {
        int_max_value = args[2].val.i;
      }
      break;
    }
    }

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    histo = ihisto->value_select(error , int_min_value , int_max_value , keep);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::FREQUENCY_DISTRIBUTION , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ValueSelect");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    int nb_variable , variable , offset;
    double real_min_value , real_max_value;
    const Vectors *ivec;
    const Sequences *iseq;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    }

    if (nb_variable == 1) {
      offset = 1;
      variable = 1;
    }

    else {
      offset = 2;

      if (args[1].tag() != AMObjType::INTEGER) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , 2 ,
                    args[1].tag.string().data() , "INT");
      }
      else {
        variable = args[1].val.i;
      }
    }

    CHECKCONDVA((nb_required - offset == 1) || (nb_required - offset == 2) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

    if (args[offset].tag() == AMObjType::INTEGER) {
      int_min_value = args[offset].val.i;

      switch (nb_required - offset) {

      case 1 : {
        int_max_value = int_min_value;
        break;
      }

      case 2 : {
        if (args[offset + 1].tag() != AMObjType::INTEGER) {
          status = false;
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , offset + 2 ,
                      args[offset + 1].tag.string().data() , "INT");
        }
        else {
          int_max_value = args[offset + 1].val.i;
        }
        break;
      }
      }
    }

    else if (args[offset].tag() == AMObjType::REAL) {
      CHECKCONDVA(nb_required - offset == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

      real_min_value = args[offset].val.r;

      if (args[offset + 1].tag() != AMObjType::REAL) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , offset + 2 ,
                    args[offset + 1].tag.string().data() , "REAL");
      }
      else {
        real_max_value = args[offset + 1].val.r;
      }
    }

    else {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ValueSelect" , offset + 1 ,
                  args[offset].tag.string().data() , "INT or REAL");
    }

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
      Vectors *vec;


      if (args[offset].tag() == AMObjType::INTEGER) {
        vec = ivec->value_select(error , &AMLOUTPUT , variable ,
                                 int_min_value , int_max_value , keep);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        vec = ivec->value_select(error , &AMLOUTPUT , variable ,
                                 real_min_value , real_max_value , keep);
      }

      if (vec) {
        STAT_model* model = new STAT_model(vec);
        return AMObj(AMObjType::VECTORS , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ValueSelect");
        return AMObj(AMObjType::ERROR);
      }
    }

    if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
        (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
        (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      Sequences *seq;
      MarkovianSequences *markovian_seq;


      if (args[offset].tag() == AMObjType::INTEGER) {
        seq = iseq->value_select(error , &AMLOUTPUT , variable ,
                                 int_min_value , int_max_value , keep);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        seq = iseq->value_select(error , &AMLOUTPUT , variable ,
                                 real_min_value , real_max_value , keep);
      }

      if (seq) {
        markovian_seq = seq->markovian_sequences(error);
        if (markovian_seq) {
          delete seq;
          STAT_model* model = new STAT_model(markovian_seq);
          return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
        }
        else {
          AMLOUTPUT << "\n";
          error.ascii_write(AMLOUTPUT , WARNING);
          STAT_model* model = new STAT_model(seq);
          return AMObj(AMObjType::SEQUENCES , model);
        }
      }

      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ValueSelect");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

//  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "ValueSelect" , 1 ,
  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_sds) , "ValueSelect" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Changement d'unite d'une variable.
 *
 *--------------------------------------------------------------*/

AMObj STAT_VariableScaling(const AMObjVector &args)

{
  bool status = true;
  int nb_variable , variable , offset;
  const Vectors *ivec;
  const Sequences *iseq;
  StatError error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "VariableScaling"));

  switch (args[0].tag()) {
  case AMObjType::VECTORS :
    ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = ivec->get_nb_variable();
    break;
  case AMObjType::MIXTURE_DATA :
    ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = ivec->get_nb_variable();
    break;
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , 1 , args[0].tag.string().data() ,
                "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  CHECKCONDVA(args.length() == offset + 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "VariableScaling" , offset + 1));

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    Vectors *vec;


    if (args[offset].tag() == AMObjType::INTEGER) {
      vec = ivec->scaling(error , variable , args[offset].val.i);
    }
    else if (args[offset].tag() == AMObjType::REAL) {
      vec = ivec->scaling(error , variable , args[offset].val.r);
    }
    else {
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , offset + 1 ,
                  args[offset].tag.string().data() , "INT or REAL");
      return AMObj(AMObjType::ERROR);
    }

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "VariableScaling");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    if (args[offset].tag() == AMObjType::INTEGER) {
      seq = iseq->scaling(error , variable , args[offset].val.i);
    }
    else if (args[offset].tag() == AMObjType::REAL) {
      seq = iseq->scaling(error , variable , args[offset].val.r);
    }
    else {
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , offset + 1 ,
                  args[offset].tag.string().data() , "INT or REAL");
      return AMObj(AMObjType::ERROR);
    }

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "VariableScaling");
      return AMObj(AMObjType::ERROR);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Arrondi des valeurs d'une variable reelle.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Round(const AMObjVector &args)

{
  string *pstr;
  bool status = true , variable_option = false , mode_option = false;
  int i;
  int nb_required , variable = I_DEFAULT;
  rounding mode = ROUND;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Round"));

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Round" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Variable") {
        if (variable_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Round" , nb_required + i + 1);
        }

        else {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Round" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else if (*pstr == "Mode") {
        if (mode_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Round" , nb_required + i + 1);
        }

        else {
          mode_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Round" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "Floor") {
              mode = FLOOR;
            }
            else if (*pstr == "Round") {
              mode = ROUND;
            }
            else if (*pstr == "Ceil") {
              mode = CEIL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(ROUND_MODE_sds) , "Round" , nb_required + i + 1 ,
                          "Floor or Round or Ceil");
            }
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Round" , nb_required + i + 1 ,
                    "Variable or Mode");
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors *ivec;
    Vectors *vec;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    vec = ivec->round(error , variable , mode);

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Round");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences *iseq;
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = iseq->round(error , variable , mode);

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Round");
      return AMObj(AMObjType::ERROR);
    }
  }

  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Round" , args[0].tag.string().data() ,
              "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Log-transform of values of vectors or sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_LogTransform(const AMObjVector &args)

{
  string *pstr;
  bool status = true , variable_option = false , base_option = false;
  int i;
  int nb_required , variable = I_DEFAULT;
  log_base base = NATURAL;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "LogTransform"));

  // argument optionnel

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LogTransform" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Variable") {
        if (variable_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "LogTransform" , nb_required + i + 1);
        }

        else {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LogTransform" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else if (*pstr == "Base") {
        if (base_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "LogTransform" , nb_required + i + 1);
        }

        else {
          base_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LogTransform" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "Natural") {
              base = NATURAL;
            }
            else if (*pstr == "Two") {
              base = TWO;
            }
            else if (*pstr == "Ten") {
              base = TEN;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(LOG_BASE_sds) , "LogTransform" , nb_required + i + 1 ,
                          "Natural or Two or Ten");
            }
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "LogTransform" , nb_required + i + 1 ,
                    "Variable or Base");
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  // argument obligatoire

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors *ivec;
    Vectors *vec;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    vec = ivec->log_transform(error , variable , base);

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "LogTransform");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences *iseq;
    Sequences *seq;
    MarkovianSequences *markovian_seq;

    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    default :
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "LogTransform" , args[0].tag.string().data() ,
                  "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
      break;
    }

    seq = iseq->log_transform(error , variable , base);

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "LogTransform");
      return AMObj(AMObjType::ERROR);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Changement du pas de regroupement de l'histogramme marginal et
 *  des histogrammes d'observation pour une variable donnee.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SelectBinWidth(const AMObjVector &args)

{
  bool status = true;
  int nb_variable , variable , offset;
  double step , min_value = D_INF;
  Vectors *vec;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SelectBinWidth"));

  switch (args[0].tag()) {
  case AMObjType::VECTORS :
    vec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = vec->get_nb_variable();
    break;
  case AMObjType::MIXTURE_DATA :
    vec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = vec->get_nb_variable();
    break;
  case AMObjType::SEQUENCES :
    seq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = seq->get_nb_variable();
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    markovian_seq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = markovian_seq->get_nb_variable();
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    markovian_seq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = markovian_seq->get_nb_variable();
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    markovian_seq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = markovian_seq->get_nb_variable();
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    markovian_seq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = markovian_seq->get_nb_variable();
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectBinWidth" , 1 , args[0].tag.string().data() ,
                "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectBinWidth" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  CHECKCONDVA((args.length() == offset + 1) || (args.length() == offset + 3) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SelectBinWidth"));

  // argument obligatoire

  switch (args[offset].tag()) {
  case AMObjType::INTEGER :
    step = args[offset].val.i;
    break;
  case AMObjType::REAL :
    step = args[offset].val.r;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectBinWidth" , offset + 1 ,
                args[offset].tag.string().data() , "INT or REAL");
  }

  // argument optionnel

  if (args.length() == offset + 3) {
    if (args[offset + 1].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectBinWidth" , offset + 2 ,
                  args[offset + 1].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[offset + 1].val.p) != "MinValue") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SelectBinWidth" , offset + 2 ,
                    "MinValue");
      }
    }

    switch (args[offset + 2].tag()) {
    case AMObjType::INTEGER :
      min_value = args[offset + 2].val.i;
      break;
    case AMObjType::REAL :
      min_value = args[offset + 2].val.r;
      break;
    default :
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectBinWidth" , offset + 2 ,
                  args[offset + 2].tag.string().data() , "INT or REAL");
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    status = vec->select_bin_width(error , variable , step , min_value);

    if (status) {
      return AMObj(AMObjType::VOID);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectBinWidth");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::SEQUENCES) {
    status = seq->select_bin_width(error , variable , step , min_value);

    if (status) {
      return AMObj(AMObjType::VOID);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectBinWidth");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    status = markovian_seq->select_bin_width(error , variable , step , min_value);

    if (status) {
      return AMObj(AMObjType::VOID);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectBinWidth");
      return AMObj(AMObjType::ERROR);
    }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection par l'identificateur.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SelectIndividual(const AMObjVector &args)

{
  string *pstr;
  bool status = true , keep = true;
  int i;
  int nb_required , nb_pattern = I_DEFAULT , *identifier = NULL;
  vector<int> vec_identifier;
  StatError error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SelectIndividual"));

  // argument obligatoire

  if (args[1].tag() != AMObjType::ARRAY) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , 2 ,
                args[1].tag.string().data() , "ARRAY");
  }
  else {
    identifier = buildIntArray(args , 1 , "SelectIndividual" , 2 , nb_pattern);
    if (!identifier) {
      status = false;
    }
    else {
      for (i = 0;i < nb_pattern;i++) {
        vec_identifier.push_back(identifier[i]);
      }
    }
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SelectIndividual" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[nb_required + 1].val.p;
      if (*pstr == "Keep") {
        keep = true;
      }
      else if (*pstr == "Reject") {
        keep = false;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(ARG_sds) , "SelectIndividual" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if (!status) {
    delete [] identifier;
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors *ivec;
    Vectors *vec;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    vec = ivec->select_individual(error , nb_pattern , identifier , keep);
    vec = ivec->select_individual(error , nb_pattern , vec_identifier , keep);
    delete [] identifier;

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectIndividual");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences *iseq;
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    seq = iseq->select_individual(error , nb_pattern , identifier , keep);
    seq = iseq->select_individual(error , nb_pattern , vec_identifier , keep);
    delete [] identifier;

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectIndividual");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::DISTANCE_MATRIX) {
    DistanceMatrix *dist_matrix;


//    dist_matrix = ((DistanceMatrix*)((STAT_model*)args[0].val.p)->pt)->select_individual(error , nb_pattern ,
//                                                                                         identifier , keep);
    dist_matrix = ((DistanceMatrix*)((STAT_model*)args[0].val.p)->pt)->select_individual(error , nb_pattern ,
                                                                                         vec_identifier , keep);
    delete [] identifier;

    if (dist_matrix) {
      STAT_model* model = new STAT_model(dist_matrix);
      return AMObj(AMObjType::DISTANCE_MATRIX , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectIndividual");
      return AMObj(AMObjType::ERROR);
    }
  }

  delete [] identifier;
//   genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , 1 , args[0].tag.string().data() ,
//               "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA or TOPS or DISTANCE_MATRIX");
   genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , 1 , args[0].tag.string().data() ,
              "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA or DISTANCE_MATRIX");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SelectVariable(const AMObjVector &args)

{
  string *pstr;
  bool status = true , keep = true;
  int i;
  int nb_required , nb_variable = I_DEFAULT , *variable = NULL;
  vector<int> vec_variable;
  StatError error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SelectVariable"));

  // argument obligatoire

  switch (args[1].tag()) {

  case AMObjType::INTEGER : {
    nb_variable = 1;
    variable = new int[nb_variable];
    variable[0] = args[1].val.i;
    vec_variable.push_back(variable[0]);
    break;
  }

  case AMObjType::ARRAY : {
    variable = buildIntArray(args , 1 , "SelectVariable" , 2 , nb_variable);
    if (!variable) {
      status = false;
    }
    else {
      for (i = 0;i < nb_variable;i++) {
        vec_variable.push_back(variable[i]);
      }
    }
    break;
  }

  default : {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectVariable" , 2 ,
                args[1].tag.string().data() , "INT or ARRAY");
    break;
  }
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectVariable" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SelectVariable" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectVariable" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[nb_required + 1].val.p;
      if (*pstr == "Keep") {
        keep = true;
      }
      else if (*pstr == "Reject") {
        keep = false;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(ARG_sds) , "SelectVariable" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if (!status) {
    delete [] variable;
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors *ivec;
    Vectors *vec;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    vec = ivec->select_variable(error , nb_variable , variable , keep);
    vec = ivec->select_variable(error , nb_variable , vec_variable , keep);
    delete [] variable;

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::SEQUENCES) {
    Sequences *seq;
    MarkovianSequences *markovian_seq;


//    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->select_variable(error , nb_variable ,
//                                                                          variable , keep);
    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->select_variable(error , nb_variable ,
                                                                          vec_variable , keep);
    delete [] variable;

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const MarkovianSequences *iseq;
    MarkovianSequences *seq;


    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    seq = iseq->select_variable(error , nb_variable , variable , keep);
    seq = iseq->select_variable(error , nb_variable , vec_variable , keep);
    delete [] variable;

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  delete [] variable;
  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectVariable" , 1 , args[0].tag.string().data() ,
              "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Sommation de variables.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SumVariable(const AMObjVector &args)

{
  bool status = true;
  int i;
  int nb_variable = I_DEFAULT , *variable = NULL;
  vector<int> vec_variable;
  StatError error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "SumVariable" , 2));

  // argument obligatoire

  if (args[1].tag() == AMObjType::ARRAY) {
    variable = buildIntArray(args , 1 , "SumVariable" , 2 , nb_variable);
    if (!variable) {
      status = false;
    }
    else {
      for (i = 0;i < nb_variable;i++) {
        vec_variable.push_back(variable[i]);
      }
    }
  }

  else {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SumVariable" , 2 ,
                args[1].tag.string().data() , "ARRAY");
  }

  if (!status) {
    delete [] variable;
    return AMObj(AMObjType::ERROR);
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors *ivec;
    Vectors *vec;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MIXTURE_DATA :
      ivec = (MixtureData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    vec = ivec->sum_variable(error , nb_variable , variable);
    vec = ivec->sum_variable(error , nb_variable , vec_variable);
    delete [] variable;

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SumVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences *iseq;
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

//    seq = iseq->sum_variable(error , nb_variable , variable);
    seq = iseq->sum_variable(error , nb_variable , vec_variable);
    delete [] variable;

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SumVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  delete [] variable;
  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SumVariable" , 1 , args[0].tag.string().data() ,
              "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables.
 *
 *--------------------------------------------------------------*/

AMObj STAT_MergeVariable(const AMObjVector &args)

{
  int i;
  bool status = true;
  int nb_required , ref_sample = I_DEFAULT;
  StatError error;


  nb_required = nb_required_computation(args);

  CHECKCONDVA(nb_required >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "MergeVariable"));

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "MergeVariable"));

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "RefSample") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "MergeVariable" , nb_required + 1 ,
                    "RefSample");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      ref_sample = args[nb_required + 1].val.i;
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::MIXTURE_DATA)) {
    const Vectors **pvec;
    vector<Vectors> ivec;
    Vectors *vec;
    StatError error;


    pvec = new const Vectors*[nb_required];
    for (i = 1;i < nb_required;i++) {
      pvec[i] = NULL;
    }

    for (i = 0;i < nb_required;i++) {
      switch (args[i].tag()) {
      case AMObjType::VECTORS :
        pvec[i] = (Vectors*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MIXTURE_DATA :
        pvec[i] = (MixtureData*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , i + 1 ,
                    args[i].tag.string().data() , "VECTORS or MIXTURE_DATA");
        break;
      }

      if ((i > 0) && (pvec[i])) {
        ivec.push_back(*pvec[i]);
      }
    }

    if (status) {
//      vec = pvec[0]->merge_variable(error , nb_required - 1 , pvec + 1 , ref_sample);
      vec = pvec[0]->merge_variable(error , nb_required - 1 , ivec , ref_sample);
    }

    delete [] pvec;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (vec) {
      STAT_model* model = new STAT_model(vec);
      return AMObj(AMObjType::VECTORS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "MergeVariable");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    const Sequences **pseq;
    vector<Sequences> iseq;
    Sequences *seq;
    const MarkovianSequences **pmarkovian_seq;
    vector<MarkovianSequences> imarkovian_seq;
    MarkovianSequences *markovian_seq;


    pseq = new const Sequences*[nb_required];
    for (i = 1;i < nb_required;i++) {
      pseq[i] = NULL;
    }

    pmarkovian_seq = new const MarkovianSequences*[nb_required];
    for (i = 0;i < nb_required;i++) {
      pmarkovian_seq[i] = NULL;
    }

    for (i = 0;i < nb_required;i++) {
      switch (args[i].tag()) {
      case AMObjType::SEQUENCES :
        pseq[i] = (Sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MARKOVIAN_SEQUENCES :
        pseq[i] = (MarkovianSequences*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (MarkovianSequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
        pseq[i] = (VariableOrderMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (VariableOrderMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::SEMI_MARKOV_DATA :
        pseq[i] = (SemiMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (SemiMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
        pseq[i] = (NonhomogeneousMarkovData*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (NonhomogeneousMarkovData*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , i + 1 , args[i].tag.string().data() ,
                    "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
        break;
      }

      if (i > 0) {
        if (pseq[i]) {
          iseq.push_back(*pseq[i]);
        }
        if (pmarkovian_seq[i]) {
          imarkovian_seq.push_back(*pmarkovian_seq[i]);
        }
      }
    }

    if (!status) {
      delete [] pseq;
      delete [] pmarkovian_seq;
      return AMObj(AMObjType::ERROR);
    }

    for (i = 0;i < nb_required;i++) {
      if (!pmarkovian_seq[i]) {
        break;
      }
    }

    if (i < nb_required) {
//      seq = pseq[0]->merge_variable(error , nb_required - 1 , pseq + 1 , ref_sample);
      seq = pseq[0]->merge_variable(error , nb_required - 1 , iseq , ref_sample);
      delete [] pseq;
      delete [] pmarkovian_seq;

      if (seq) {
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "MergeVariable");
        return AMObj(AMObjType::ERROR);
      }
    }

    else {
//      markovian_seq = pmarkovian_seq[0]->merge_variable(error , nb_required - 1 ,
//                                                        pmarkovian_seq + 1 , ref_sample);
      markovian_seq = pmarkovian_seq[0]->merge_variable(error , nb_required - 1 ,
                                                        imarkovian_seq , ref_sample);
      delete [] pseq;
      delete [] pmarkovian_seq;

      if (markovian_seq) {
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "MergeVariable");
        return AMObj(AMObjType::ERROR);
      }
    }
  }

  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , 1 , args[0].tag.string().data() ,
              "VECTORS or MIXTURE_DATA or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Selection d'echantillons {temps, nombre d'evenements}
 *  sur un critere de nombre d'evenements.
 *
 *--------------------------------------------------------------*/

AMObj STAT_NbEventSelect(const AMObjVector &args)

{
  bool status = true;
  const TimeEvents *itimev;
  TimeEvents *timev;
  StatError error;


  CHECKCONDVA(args.length() == 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "NbEventSelect" , 3));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (TimeEvents*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (RenewalData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "NbEventSelect" , 1 ,
                args[0].tag.string().data() , "TIME_EVENTS or RENEWAL_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "NbEventSelect" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (args[2].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "NbEventSelect" , 3 ,
                args[2].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  timev = itimev->nb_event_select(error , args[1].val.i , args[2].val.i);

  if (timev) {
    STAT_model* model = new STAT_model(timev);
    return AMObj(AMObjType::TIME_EVENTS , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "NbEventSelect");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Changement de l'unite de temps d'une structure de donnees renouvellement.
 *
 *--------------------------------------------------------------*/

AMObj STAT_TimeScaling(const AMObjVector &args)

{
  bool status = true;
  const TimeEvents *itimev;
  TimeEvents *timev;
  StatError error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "TimeScaling" , 2));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (TimeEvents*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (RenewalData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TimeScaling" , 1 ,
                args[0].tag.string().data() , "TIME_EVENTS or RENEWAL_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TimeScaling" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  timev = itimev->time_scaling(error , args[1].val.i);

  if (timev) {
    STAT_model* model = new STAT_model(timev);
    return AMObj(AMObjType::TIME_EVENTS , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "TimeScaling");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection d'echantillons {temps, nombre d'evenements}
 *  sur un critere de temps d'observation.
 *
 *--------------------------------------------------------------*/

AMObj STAT_TimeSelect(const AMObjVector &args)

{
  bool status = true;
  int min_time , max_time;
  const TimeEvents *itimev;
  TimeEvents *timev;
  StatError error;


  CHECKCONDVA((args.length() == 2) || (args.length() == 3) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "TimeSelect"));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (TimeEvents*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (RenewalData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TimeSelect" , 1 ,
                args[0].tag.string().data() , "TIME_EVENTS or RENEWAL_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TimeSelect" , 2 ,
                args[1].tag.string().data() , "INT");
  }
  else {
    min_time = args[1].val.i;
  }

  switch (args.length()) {

  case 2 : {
    max_time = min_time;
    break;
  }

  case 3 : {
    if (args[2].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TimeSelect" , 3 ,
                  args[2].tag.string().data() , "INT");
    }
    else {
      max_time = args[2].val.i;
    }
    break;
  }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  timev = itimev->time_select(error , min_time , max_time);

  if (timev) {
    STAT_model* model = new STAT_model(timev);
    return AMObj(AMObjType::TIME_EVENTS , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "TimeSelect");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection des sequences correspondant a une plage de longueurs.
 *
 *--------------------------------------------------------------*/

AMObj STAT_LengthSelect(const AMObjVector &args)

{
  string *pstr;
  bool status = true , keep = true;
  int nb_required , min_length , max_length;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = nb_required_computation(args);

  CHECKCONDVA((nb_required == 2) || (nb_required == 3) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "LengthSelect"));

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "LengthSelect"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LengthSelect" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LengthSelect" , 2 ,
                args[1].tag.string().data() , "INT");
  }
  else {
    min_length = args[1].val.i;
  }

  switch (nb_required) {

  case 2 : {
    max_length = min_length;
    break;
  }

  case 3 : {
    if (args[2].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LengthSelect" , 3 ,
                  args[2].tag.string().data() , "INT");
    }
    else {
      max_length = args[2].val.i;
    }
    break;
  }
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LengthSelect" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "LengthSelect" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "LengthSelect" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[nb_required + 1].val.p;
      if (*pstr == "Keep") {
        keep = true;
      }
      else if (*pstr == "Reject") {
        keep = false;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(ARG_sds) , "LengthSelect" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->length_select(error , &AMLOUTPUT , min_length , max_length , keep);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "LengthSelect");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Differences entre donnees et residus.
 *
 *--------------------------------------------------------------*/

AMObj STAT_DifferenceVariable(const AMObjVector &args)

{
  bool status = true;
  const Sequences *iseq , *residual;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "DifferenceVariable" , 2));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "DifferenceVariable" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  switch (args[1].tag()) {
  case AMObjType::SEQUENCES :
    residual = (Sequences*)((STAT_model*)args[1].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    residual = (MarkovianSequences*)((STAT_model*)args[1].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "DifferenceVariable" , 2 , args[1].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES");
    break;
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->difference_variable(error , *residual);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "DifferenceVariable");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Decalage d'une variable.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ShiftVariable(const AMObjVector &args)

{
  bool status = true;
  int nb_variable , variable , offset , lag;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ShiftVariable"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ShiftVariable" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  nb_variable = iseq->get_nb_variable();

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ShiftVariable" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  CHECKCONDVA(args.length() == offset + 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ShiftVariable", offset + 1));

  if (args[offset].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ShiftVariable" , offset + 1 ,
                args[offset].tag.string().data() , "INT");
  }
  else {
    lag = args[offset].val.i;
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->shift_variable(error , variable , lag);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ShiftVariable");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Suppression des premieres/dernieres series d'une valeur donne.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RemoveRun(const AMObjVector &args)

{
  string *pstr;
  bool status = true;
  int nb_required , nb_variable , variable , offset , value , max_run_length = I_DEFAULT;
  run_position position;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() >= 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RemoveRun"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_variable = iseq->get_nb_variable();

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  nb_required = offset + 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RemoveRun"));

  if (args[offset].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , offset + 1 ,
                args[offset].tag.string().data() , "INT");
  }
  else {
    value = args[offset].val.i;
  }

  if (args[offset + 1].tag() != AMObjType::STRING) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , offset + 2 ,
                args[offset + 1].tag.string().data() , "STRING");
  }
  else {
    pstr = (AMString*)args[offset + 1].val.p;
    if (*pstr == "Begin") {
      position = BEGIN_RUN;
    }
    else if (*pstr == "End") {
      position = END_RUN;
    }
    else {
      status = false;
      genAMLError(ERRORMSG(ARG_sds) , "RemoveRun" , offset + 2 , "Begin or End");
    }
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "MaxLength") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "RemoveRun" , nb_required + 1 ,
                    "MaxLength");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveRun" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      max_run_length = args[nb_required + 1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->remove_run(error , variable , value , position , max_run_length);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveRun");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Inversion du sens de parcours de sequences ou de cimes.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Reverse(const AMObjVector &args)

{
  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Reverse"));

  if (args[0].tag() == AMObjType::SEQUENCES) {
    Sequences *seq;
    StatError error;


    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->reverse(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Reverse");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    MarkovianSequences *iseq , *seq;


    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = new MarkovianSequences(*iseq , REVERSE);

    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }

//  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Reverse" , args[0].tag.string().data() ,
//              "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA or TOPS");
  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Reverse" , args[0].tag.string().data() ,
              "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RemoveIndexParameter(const AMObjVector &args)

{
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "RemoveIndexParameter"));

  // argument obligatoire

  switch (args[0].tag()) {

  case AMObjType::SEQUENCES : {
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

    if (seq) {
      markovian_seq = seq->markovian_sequences(error);
      if (markovian_seq) {
        delete seq;
        STAT_model* model = new STAT_model(markovian_seq);
        return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
      }
      else {
        AMLOUTPUT << "\n";
        error.ascii_write(AMLOUTPUT , WARNING);
        STAT_model* model = new STAT_model(seq);
        return AMObj(AMObjType::SEQUENCES , model);
      }
    }

    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::MARKOVIAN_SEQUENCES : {
    MarkovianSequences *seq;


    seq = ((MarkovianSequences*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::VARIABLE_ORDER_MARKOV_DATA : {
    VariableOrderMarkovData *seq;


    seq = ((VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::VARIABLE_ORDER_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::SEMI_MARKOV_DATA : {
    SemiMarkovData *seq;


    seq = ((SemiMarkovData*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEMI_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA : {
    NonhomogeneousMarkovData *seq;


    seq = ((NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::NONHOMOGENEOUS_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  default : {
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "RemoveIndexParameter" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie d'un objet Sequences avec transformation du parametre d'index implicite
 *  en parametre d'index explicite.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExplicitIndexParameter(const AMObjVector &args)

{
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ExplicitIndexParameter"));

  // argument obligatoire

  switch (args[0].tag()) {

  case AMObjType::SEQUENCES : {
    Sequences *seq;
    MarkovianSequences *markovian_seq;


    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->explicit_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExplicitIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::MARKOVIAN_SEQUENCES : {
    MarkovianSequences *seq;


    seq = ((MarkovianSequences*)((STAT_model*)args[0].val.p)->pt)->explicit_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExplicitIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::VARIABLE_ORDER_MARKOV_DATA : {
    VariableOrderMarkovData *seq;


    seq = ((VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt)->explicit_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::VARIABLE_ORDER_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExplicitIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::SEMI_MARKOV_DATA : {
    SemiMarkovData *seq;


    seq = ((SemiMarkovData*)((STAT_model*)args[0].val.p)->pt)->explicit_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEMI_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExplicitIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA : {
    NonhomogeneousMarkovData *seq;


    seq = ((NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt)->explicit_index_parameter(error);

    if (seq) {
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::NONHOMOGENEOUS_MARKOV_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExplicitIndexParameter");
      return AMObj(AMObjType::ERROR);
    }
  }

  default : {
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ExplicitIndexParameter" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }
  }
}


/*--------------------------------------------------------------*
 *
 *  Selection sur les valeurs prises par le parametre d'index.
 *
 *--------------------------------------------------------------*/

AMObj STAT_IndexParameterSelect(const AMObjVector &args)

{
  string *pstr;
  bool status = true , keep = true;
  int nb_required , min_index_parameter , max_index_parameter;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = nb_required_computation(args);

  CHECKCONDVA((nb_required == 2) || (nb_required == 3) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "IndexParameterSelect"));

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "IndexParameterSelect"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterSelect" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterSelect" , 2 ,
                args[1].tag.string().data() , "INT");
  }
  else {
    min_index_parameter = args[1].val.i;
  }

  switch (nb_required) {

  case 2 : {
    max_index_parameter = min_index_parameter;
    break;
  }

  case 3 : {
    if (args[2].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterSelect" , 3 ,
                  args[2].tag.string().data() , "INT");
    }
    else {
      max_index_parameter = args[2].val.i;
    }
    break;
  }
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterSelect" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "IndexParameterSelect" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterSelect" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "STRING");
    }
    else {
      pstr = (AMString*)args[nb_required + 1].val.p;
      if (*pstr == "Keep") {
        keep = true;
      }
      else if (*pstr == "Reject") {
        keep = false;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(ARG_sds) , "IndexParameterSelect" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->index_parameter_select(error , &AMLOUTPUT , min_index_parameter ,
                                     max_index_parameter , keep);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "IndexParameterSelect");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Troncature de sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Truncate(const AMObjVector &args)

{
  bool status = true;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() == 2 , genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Truncate", 2));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Truncate" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Truncate" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->truncate(error , args[1].val.i);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Truncate");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction de sous-sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_IndexParameterExtract(const AMObjVector &args)

{
  bool status = true;
  int nb_required , max_index_parameter = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "IndexParameterExtract"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterExtract" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterExtract" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterExtract" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "MaxIndex") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "IndexParameterExtract" , nb_required + 1 ,
                    "MaxIndex");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "IndexParameterExtract" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      max_index_parameter = args[nb_required + 1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->index_parameter_extract(error , args[1].val.i , max_index_parameter);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "IndexParameterExtract");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Extraction par segmentation.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SegmentationExtract(const AMObjVector &args)

{
  string *pstr;
  bool status = true , mode_option = false , keep = true , concatenation_option = false ,
       concatenation = false;
  int i;
  int nb_required , nb_value = I_DEFAULT , *value = NULL;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 3;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SegmentationExtract"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  switch (args[2].tag()) {

  case AMObjType::INTEGER : {
    nb_value = 1;
    value = new int[nb_value];
    value[0] = args[2].val.i;
    break;
  }

  case AMObjType::ARRAY :{
    value = buildIntArray(args , 2 , "SegmentationExtract" , 3 , nb_value);
    if (!value) {
      status = false;
    }
    break;
  }

  default : {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , 3 ,
                args[2].tag.string().data() , "INT or ARRAY");
    break;
  }
  }

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Mode") {
        if (mode_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "SegmentationExtract" , nb_required + i + 1);
        }

        else {
          mode_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }

          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "Keep") {
              keep = true;
            }
            else if (*pstr == "Reject") {
              keep = false;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(ARG_sds) , "SegmentationExtract" ,
                          nb_required + i + 1 , "Keep or Reject");
            }
          }
        }
      }

      else if (*pstr == "Concatenation") {
        if (concatenation_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "SegmentationExtract" , nb_required + i + 1);
        }

        else {
          concatenation_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            concatenation = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SegmentationExtract" , nb_required + i + 1 ,
                    "Mode or Concatenation");
      }
    }
  }

  if (status) {
    seq = iseq->segmentation_extract(error , args[1].val.i , nb_value , value ,
                                     keep , concatenation);
  }

  delete [] value;

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "SegmentationExtract");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Cumul des valeurs successives des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Cumulate(const AMObjVector &args)

{
  bool status = true;
  int nb_required , variable = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cumulate"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Cumulate" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cumulate" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Variable") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Cumulate" , nb_required + 1 , "Variable");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cumulate" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      variable = args[nb_required + 1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->cumulate(error , variable);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Cumulate");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Differenciation au 1er ordre des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Difference(const AMObjVector &args)

{
  string *pstr;
  bool status = true , variable_option = false , first_element_option = false , first_element = false;
  int i;
  int nb_required , variable = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Difference"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Difference" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Difference" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Variable") {
        if (variable_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Difference" , nb_required + i + 1);
        }

        else {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Difference" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else if (*pstr == "FirstElement") {
        if (first_element_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Difference" , nb_required + i + 1);
        }

        else {
          first_element_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Difference" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            first_element = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Difference" , nb_required + i + 1 ,
                    "Variable or FirstElement");
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->difference(error , variable , first_element);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Difference");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Computation of relative growth rates on the basis of cumulative dimensions.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RelativeGrowthRate(const AMObjVector &args)

{
  string *pstr;
  bool status = true;
  int nb_required;
  double growth_factor = GROWTH_FACTOR;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RelativeGrowthRate"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "RelativeGrowthRate" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RelativeGrowthRate" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "GrowthFactor") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "RelativeGrowthRate" , nb_required + 1 , "GrowthFactor");
      }
    }

    if (args[nb_required + 1].tag() == AMObjType::INTEGER) {
      growth_factor = args[nb_required + 1].val.i;
    }
    if (args[nb_required + 1].tag() == AMObjType::REAL) {
      growth_factor = args[nb_required + 1].val.r;
    }
    else {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RelativeGrowthRate" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT or REAL");
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->relative_growth_rate(error , growth_factor);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "RelativeGrowthRate");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Normalisation des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SequenceNormalization(const AMObjVector &args)

{
  bool status = true;
  int nb_required , variable = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SequenceNormalization"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "SequenceNormalization" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SequenceNormalization" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Variable") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SequenceNormalization" , nb_required + 1 , "Variable");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SequenceNormalization" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      variable = args[nb_required + 1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->sequence_normalization(error , variable);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "SequenceNormalization");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Filtrage de type moyenne mobile des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_MovingAverage(const AMObjVector &args)

{
  string *pstr;
  bool status = true , variable_option = false , begin_end_option = false , begin_end = false ,
       segmentation_option = false , segmentation = false , output_option = false;
  int i;
  int nb_required , nb_point = I_DEFAULT , int_sum , variable = I_DEFAULT , *int_filter;
  sequence_type output = TREND;
  double sum , *filter;
  vector<double> vec_filter;
  const Distribution *dist;
  const Sequences *iseq = NULL;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) || (args.length() == nb_required + 6) ||
              (args.length() == nb_required + 8) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "MovingAverage"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "MovingAverage" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  filter = NULL;
  dist = NULL;

  if (args[1].tag() == AMObjType::ARRAY) {
    Array* parray = (Array*)args[1].val.p;

    if (parray->surfaceType() == AMObjType::INTEGER) {
      int_filter = buildIntArray(args , 1 , "MovingAverage" , 2 , nb_point);
      if (!int_filter) {
        status = false;
      }

      else {
        nb_point--;

        int_sum = 0;
        for (i = 0;i < nb_point;i++) {
          int_sum += 2 * int_filter[i];
        }
        int_sum += int_filter[i];

        if (int_sum == 0) {
          status = false;
          genAMLError(ERRORMSG(FILTER_NORMALIZATION_sd) , "MovingAverage" , 2);
        }

        else {
          filter = new double[2 * nb_point + 1];
          for (i = 0;i < nb_point;i++) {
            filter[i] = (double)int_filter[i] / (double)int_sum;
            filter[2 * nb_point - i] = filter[i];
          }
          filter[i] = (double)int_filter[i] / (double)int_sum;
        }

        delete [] int_filter;
      }
    }

    else if (parray->surfaceType() == AMObjType::REAL) {
      filter = buildRealArray(args , 1 , "MovingAverage" , 2 , nb_point , true);
      if (!filter) {
        status = false;
      }

      else {
        sum = 0.;
        for (i = 0;i < nb_point;i++) {
          sum += 2 * filter[i];
          filter[2 * nb_point - i] = filter[i];
        }
        sum += filter[i];

        if ((sum < 1. - STAT_FILTER_ROUNDNESS) || (sum > 1. + STAT_FILTER_ROUNDNESS)) {
          status = false;
          genAMLError(ERRORMSG(FILTER_NORMALIZATION_sd) , "MovingAverage" , 2);
        }
      }
    }

    else {
      status = false;
      genAMLError(ERRORMSG(ARRAY_ELEMENT_TYPE_1_sdss) , "MovingAverage" , 2 ,
                  (parray->surfaceType()).string().data() , "INT or REAL");
    }

    if (filter) {
      for (i = 0;i < 2 * nb_point + 1;i++) {
        vec_filter.push_back(filter[i]);
      }
    }
  }

  else {
    switch (args[1].tag()) {
    case AMObjType::DISTRIBUTION :
      dist = new Distribution(*((Distribution*)((DiscreteParametricModel*)((STAT_model*)args[1].val.p)->pt)));
      break;
    case AMObjType::DISCRETE_MIXTURE :
      dist = new Distribution(*((Distribution*)((DiscreteMixture*)((STAT_model*)args[1].val.p)->pt)));
      break;
    case AMObjType::CONVOLUTION :
      dist = new Distribution(*((Distribution*)((Convolution*)((STAT_model*)args[1].val.p)->pt)));
      break;
    case AMObjType::COMPOUND :
      dist = new Distribution(*((Distribution*)((Compound*)((STAT_model*)args[1].val.p)->pt)));
      break;
    default :
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , 2 , args[1].tag.string().data() ,
                  "ARRAY or DISTRIBUTION or DISCRETE_MIXTURE or CONVOLUTION or COMPOUND");
      break;
    }
  }

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Variable") {
        if (variable_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
        }

        else {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else if (*pstr == "BeginEnd") {
        if (begin_end_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
        }

        else {
          begin_end_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            begin_end = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "Segmentation") {
        if (segmentation_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
        }

        else {
          segmentation_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            segmentation = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "Output") {
        if (output_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
        }

        else {
          output_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "Sequence") {
              output = SEQUENCE;
            }
            else if (*pstr == "Trend") {
              output = TREND;
            }
            else if ((*pstr == "SubtractionResidual") || (*pstr == "Residual")) {
              output = SUBTRACTION_RESIDUAL;
            }
            else if (*pstr == "DivisionResidual") {
              output = DIVISION_RESIDUAL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(MOVING_AVERAGE_OUTPUT_sds) , "MovingAverage" , nb_required + i + 1 ,
                          "Sequence or Trend or (Subtraction)Residual or DivisionResidual");
            }
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "MovingAverage" , nb_required + i + 1 ,
                    "Variable or BeginEnd or Segmentation or Output");
      }
    }
  }

  if ((segmentation_option) && (iseq) && (iseq->get_type(0) != STATE)) {
    status = false;
    genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "MovingAverage" , "Segmentation");
  }
  if ((!begin_end) && (segmentation_option)) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_sss) , "MovingAverage" , "BeginEnd" , "Segmentation");
  }

  if (status) {
    if (!dist) {
//      seq = iseq->moving_average(error , nb_point , filter , variable , begin_end , segmentation , output);
      seq = iseq->moving_average(error , nb_point , vec_filter , variable , begin_end , segmentation , output);
    }
    else {
      seq = iseq->moving_average(error , *dist , variable , begin_end , segmentation , output);
    }
  }

  delete [] filter;
  delete dist;

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "MovingAverage");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul de la sequence des moyennes (et des ecart-types).
 *
 *--------------------------------------------------------------*/

AMObj STAT_PointwiseAverage(const AMObjVector &args)

{
  string *pstr;
  output_format format = ASCII;
//  char *file_name = NULL;
  string file_name = "";
  bool status = true , circular_option = false , circular = false , robust_option = false ,
       robust = false , dispersion_option = false , dispersion = false , output_option = false ,
       file_name_option = false , format_option = false;
  int i;
  int nb_required;
  sequence_type output = SEQUENCE;
  const Sequences *iseq;
  Sequences *seq;
  StatError error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) || (args.length() == nb_required + 6) ||
              (args.length() == nb_required + 8) || (args.length() == nb_required + 10) ||
              (args.length() == nb_required + 12) || (args.length() == nb_required + 14) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "PointwiseAverage"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "PointwiseAverage" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "Circular") {
        if (circular_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          circular_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            circular = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "Robust") {
        if (robust_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          robust_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            robust = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "Dispersion") {
        if (dispersion_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          dispersion_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            dispersion = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "Output") {
        if (output_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          output_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "Sequence") {
              output = SEQUENCE;
            }
            else if ((*pstr == "SubtractionResidual") || (*pstr == "Residual")) {
              output = SUBTRACTION_RESIDUAL;
            }
            else if (*pstr == "StandardizedResidual") {
              output = STANDARDIZED_RESIDUAL;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(POINTWISE_AVERAGE_OUTPUT_sds) , "PointwiseAverage" , nb_required + i + 1 ,
                          "Sequence or (Subtraction)Residual or StandardizedResidual");
            }
          }
        }
      }

      else if (*pstr == "FileName") {
        if (file_name_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          file_name_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            file_name = (string)((AMString*)args[nb_required + i * 2 + 1].val.p)->data();
          }
        }
      }

      else if (*pstr == "Format") {
        if (format_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
        }

        else {
          format_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "ASCII") {
              format = ASCII;
            }
            else if (*pstr == "SpreadSheet") {
              format = SPREADSHEET;
            }
            else {
              status = false;
              genAMLError(ERRORMSG(K_FILE_FORMAT_ERR_sds) , "PointwiseAverage" ,
                          nb_required + i + 1 , "ASCII or SpreadSheet");
            }
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "PointwiseAverage" , nb_required + i + 1 ,
                    "Circular or Robust or Dispersion or Output or FileName or Format");
      }
    }
  }

  if ((circular) && (output == STANDARDIZED_RESIDUAL)) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_sss) , "PointwiseAverage" , "Circular" , "Output");
  }

  if ((circular_option) && (robust_option)) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_sss) , "PointwiseAverage" , "Circular" , "Robust");
  }

  if ((format_option) && (!file_name_option)) {
    status = false;
    genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "PointwiseAverage" , "Format");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->pointwise_average(error , circular , robust , dispersion ,
                                output , file_name , format);

  if (seq) {
    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::SEQUENCES , model);
   }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "PointwiseAverage");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences de temps de retour pour une valeur prise
 *  par une variable entiere.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RecurrenceTimeSequences(const AMObjVector &args)

{
  bool status = true;
  int nb_variable , variable , offset;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RecurrenceTimeSequences"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RecurrenceTimeSequences" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_variable = iseq->get_nb_variable();

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RecurrenceTimeSequences" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  CHECKCONDVA(args.length() == offset + 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "RecurrenceTimeSequences" , offset + 1));

  if (args[offset].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RecurrenceTimeSequences" , offset + 1 ,
                args[offset].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->recurrence_time_sequences(error , variable , args[offset].val.i);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "RecurrenceTimeSequences");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Calcul des sequences de temps de sejour pour une variable entiere.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SojournTimeSequences(const AMObjVector &args)

{
  bool status = true;
  int nb_variable , variable;
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SojournTimeSequences"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SojournTimeSequences" , 1 , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_variable = iseq->get_nb_variable();

  if (nb_variable == 1) {
    CHECKCONDVA(args.length() == 1 ,
                genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "SojournTimeSequences"));

    variable = 1;
  }

  else {
    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "SojournTimeSequences" , 2));

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SojournTimeSequences" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->sojourn_time_sequences(error , variable);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "SojournTimeSequences");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Discretisation des positions.
 *
 *--------------------------------------------------------------*/

AMObj STAT_TransformPosition(const AMObjVector &args)

{
  bool status = true;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "TransformPosition" , 2));

  if (args[0].tag() != AMObjType::SEQUENCES) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TransformPosition" , 1 ,
                args[0].tag.string().data() , "SEQUENCES");
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "TransformPosition" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->transform_position(error , args[1].val.i);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "TransformPosition");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Croisement des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Cross(const AMObjVector &args)

{
  const Sequences *iseq;
  Sequences *seq;
  MarkovianSequences *markovian_seq;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Cross"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Cross" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->cross(error);

  if (seq) {
    markovian_seq = seq->markovian_sequences(error);
    if (markovian_seq) {
      delete seq;
      STAT_model* model = new STAT_model(markovian_seq);
      return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
    }
    else {
      AMLOUTPUT << "\n";
      error.ascii_write(AMLOUTPUT , WARNING);
      STAT_model* model = new STAT_model(seq);
      return AMObj(AMObjType::SEQUENCES , model);
    }
  }

  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Cross");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Copie des sequence avec extraction des temps de sejour initiaux.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ComputeInitialRun(const AMObjVector &args)

{
  const MarkovianSequences *iseq;
  MarkovianSequences *seq;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ComputeInitialRun"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ComputeInitialRun" , args[0].tag.string().data() ,
                "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->initial_run_computation(error);

  if (seq) {
    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ComputeInitialRun");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Ajout d'une serie de vecteurs absorbants a la fin de chaque sequence.
 *
 *--------------------------------------------------------------*/

AMObj STAT_AddAbsorbingRun(const AMObjVector &args)

{
  string *pstr;
  bool status = true , add_variable_option = false , add_variable = false ,
       run_length_option = false , sequence_length_option = false;
  int i;
  int nb_required , run_length = I_DEFAULT , sequence_length = I_DEFAULT;
  const MarkovianSequences *iseq;
  MarkovianSequences *seq;
  StatError error;


  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "AddAbsorbingRun"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "AddAbsorbingRun" , args[0].tag.string().data() ,
                "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "AddAbsorbingRun"));

  // arguments optionnels

  for (i = 0;i < (args.length() - nb_required) / 2;i++) {
    if (args[nb_required + i * 2].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + i + 1 ,
                  args[nb_required + i * 2].tag.string().data() , "OPTION");
    }

    else {
      pstr = (AMString*)args[nb_required + i * 2].val.p;

      if (*pstr == "AddVariable") {
        if (add_variable_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "AddAbsorbingRun" , nb_required + i + 1);
        }

        else {
          add_variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            add_variable = args[nb_required + i * 2 + 1].val.b;
          }
        }
      }

      else if (*pstr == "RunLength") {
        if (run_length_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "AddAbsorbingRun" , nb_required + i + 1);
        }

        else {
          run_length_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            run_length = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else if (*pstr == "SequenceLength") {
        if (sequence_length_option) {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "AddAbsorbingRun" , nb_required + i + 1);
        }

        else {
          sequence_length_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            sequence_length = args[nb_required + i * 2 + 1].val.i;
          }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "AddAbsorbingRun" , nb_required + i + 1 ,
                    "AddVariable or RunLength or SequenceLength");
      }
    }
  }

  if ((run_length_option) && (sequence_length_option)) {
    status = false;
    genAMLError(ERRORMSG(INCOMPATIBLE_OPTIONS_sss) , "AddAbsorbingRun" , "RunLength" , "SequenceLength");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->add_absorbing_run(error , run_length , sequence_length , add_variable);

  if (seq) {
    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "AddAbsorbingRun");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Suppression des valeurs non-representees d'une variable entiere.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ConsecutiveValues(const AMObjVector &args)

{
  bool status = true , add_variable = false;
  int nb_required , nb_variable , variable , offset;
  const MarkovianSequences *iseq;
  MarkovianSequences *seq;
  StatError error;


  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ConsecutiveValues"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ConsecutiveValues" , args[0].tag.string().data() ,
                "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_variable = iseq->get_nb_variable();

  if (nb_variable == 1) {
    offset = 1;
    variable = 1;
  }

  else {
    offset = 2;

    if (args[1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ConsecutiveValues" , 2 ,
                  args[1].tag.string().data() , "INT");
    }
    else {
      variable = args[1].val.i;
    }
  }

  nb_required = offset;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ConsecutiveValues"));

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ConsecutiveValues" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "AddVariable") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "ConsecutiveValues" , nb_required + 1 , "AddVariable");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::BOOL) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ConsecutiveValues" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "BOOL");
    }
    else {
      add_variable = args[nb_required + 1].val.b;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->consecutive_values(error , &AMLOUTPUT , variable , add_variable);

  if (seq) {
    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ConsecutiveValues");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Segmentation de sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Split(const AMObjVector &args)

{
  bool status = true;
  const MarkovianSequences *iseq;
  MarkovianSequences *seq;
  StatError error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Split" , 2));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (MarkovianSequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (VariableOrderMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (SemiMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (NonhomogeneousMarkovData*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Split" , 1 , args[0].tag.string().data() ,
                "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Split" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->split(error , args[1].val.i);

  if (seq) {
    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Split");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Conversion des resultats d'un clustering par partitionnement en
 *  matrice des distances.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToDistanceMatrix(const AMObjVector &args)

{
  DistanceMatrix *dist_matrix;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToDistanceMatrix"));

  CHECKCONDVA(args[0].tag() == AMObjType::CLUSTERS ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToDistanceMatrix" ,
                          args[0].tag.string().data() , "CLUSTERS"));

  dist_matrix = new DistanceMatrix(*((Clusters*)((STAT_model*)args[0].val.p)->pt));

  STAT_model* model = new STAT_model(dist_matrix);
  return AMObj(AMObjType::DISTANCE_MATRIX , model);
}


/*--------------------------------------------------------------*
 *
 *  Symetrisation d'une matrice des distances.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Symmetrize(const AMObjVector &args)

{
  DistanceMatrix *dist_matrix;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Symmetrize"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTANCE_MATRIX ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Symmetrize" , 1 ,
                          args[0].tag.string().data() , "DISTANCE_MATRIX"));

  dist_matrix = ((DistanceMatrix*)((STAT_model*)args[0].val.p)->pt)->symmetrize(error);

  if (dist_matrix) {
    STAT_model* model = new STAT_model(dist_matrix);
    return AMObj(AMObjType::DISTANCE_MATRIX , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Symmetrize");
    return AMObj(AMObjType::ERROR);
  }
}


/*--------------------------------------------------------------*
 *
 *  Denormalisation d'une matrice des distances.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Unnormalize(const AMObjVector &args)

{
  DistanceMatrix *dist_matrix;
  StatError error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Unnormalize"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTANCE_MATRIX ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Unnormalize" , 1 ,
                          args[0].tag.string().data() , "DISTANCE_MATRIX"));

  dist_matrix = ((DistanceMatrix*)((STAT_model*)args[0].val.p)->pt)->unnormalize(error);

  if (dist_matrix) {
    STAT_model* model = new STAT_model(dist_matrix);
    return AMObj(AMObjType::DISTANCE_MATRIX , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "Unnormalize");
    return AMObj(AMObjType::ERROR);
  }
}
