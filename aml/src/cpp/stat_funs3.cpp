/* -*-c++-*-
 *  ----------------------------------------------------------------------------
 *
 *       AMAPmod: Exploring and Modeling Plant Architecture
 *
 *       Copyright 1995-2002 UMR Cirad/Inra Modelisation des Plantes
 *
 *       File author(s): Y. Guedon (yann.guedon@cirad.fr)
 *
 *       $Source$
 *       $Id$
 *
 *       Forum for AMAPmod developers: amldevlp@cirad.fr
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

#include "stat_tool/stat_tools.h"
#include "stat_tool/distribution.h"
#include "stat_tool/mixture.h"
#include "stat_tool/convolution.h"
#include "stat_tool/compound.h"
#include "stat_tool/vectors.h"
#include "stat_tool/regression.h"
#include "stat_tool/curves.h"
#include "stat_tool/markovian.h"
#include "stat_tool/distance_matrix.h"

#include "sequence_analysis/renewal.h"
#include "sequence_analysis/sequences.h"
#include "sequence_analysis/variable_order_markov.h"
#include "sequence_analysis/hidden_variable_order_markov.h"
#include "sequence_analysis/semi_markov.h"
#include "sequence_analysis/hidden_semi_markov.h"
#include "sequence_analysis/nonhomogeneous_markov.h"
#include "sequence_analysis/tops.h"

#include "aml/ammodel.h"
#include "aml/parseraml.h"
#include "aml/kernel_err.h"
#include "aml/amstring.h"
#include "aml/array.h"

#include "aml/stat_module.h"


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
 *  Conversion d'un histogramme en loi.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToDistribution(const AMObjVector &args)

{
  Parametric_model *dist;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToDistribution"));

  CHECKCONDVA(args[0].tag() == AMObjType::HISTOGRAM ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToDistribution" ,
                          args[0].tag.string().data() , "HISTOGRAM"));

  dist = ((Distribution_data*)((STAT_model*)args[0].val.p)->pt)->extract_model(error);

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
 *  Conversion d'une loi en histogramme.
 *
 *--------------------------------------------------------------*/

AMObj STAT_ToHistogram(const AMObjVector &args)

{
  Distribution_data *histo;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToHistogram"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTRIBUTION ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToHistogram" ,
                          args[0].tag.string().data() , "DISTRIBUTION"));

  histo = ((Parametric_model*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

  if (histo) {
    STAT_model* model = new STAT_model(histo);
    return AMObj(AMObjType::HISTOGRAM , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "ToHistogram");
    return AMObj(AMObjType::ERROR);
  }
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

  if (args[0].tag() == AMObjType::MIXTURE) {
    RWCString *pstr;
    Parametric_model *dist;
    const Mixture *mixture;
    const Mixture_data *mixture_data;
    Format_error error;


    mixture = (Mixture*)((STAT_model*)args[0].val.p)->pt;
    mixture_data = mixture->get_mixture_data();

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Weight") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

      dist = new Parametric_model(*(mixture->get_weight()) ,
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

      dist = new Parametric_model(*((Distribution*)mixture) , (Histogram*)mixture_data);
//      dist = new Parametric_model((Distribution&)*mixture , (Histogram*)mixture_data);  segmentation fault

      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                "Weight or Component or Mixture");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::CONVOLUTION) {
    RWCString *pstr;
    Parametric_model *dist;
    const Convolution *convolution;
    const Convolution_data *convolution_data;
    Format_error error;


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

      dist = new Parametric_model(*((Distribution*)convolution) , (Histogram*)convolution_data);

      STAT_model* model = new STAT_model(dist);
      return AMObj(AMObjType::DISTRIBUTION , model);
    }

    genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                "Elementary or Convolution");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::COMPOUND) {
    RWCString *pstr;
    Parametric_model *dist;
    const Compound *compound;
    const Compound_data *compound_data;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

    compound = (Compound*)((STAT_model*)args[0].val.p)->pt;
    compound_data = compound->get_compound_data();

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Sum") {
      dist = new Parametric_model(*(compound->get_sum_distribution()) ,
                                  (compound_data ? compound_data->get_sum_histogram() : 0));
    }
    else if (*pstr == "Elementary") {
      dist = new Parametric_model(*(compound->get_distribution()) ,
                                  (compound_data ? compound_data->get_histogram() : 0));
    }
    else if (*pstr == "Compound") {
      dist = new Parametric_model(*((Distribution*)compound) , (Histogram*)compound_data);
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
    RWCString *pstr;
    int ident , time;
    Parametric_model *dist;
    Format_error error;


    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "NbEvent") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 3));

      ident = NB_EVENT;

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
        ident = INTER_EVENT;
      }
      else if (*pstr == "LengthBias") {
        ident = LENGTH_BIAS;
      }
      else if (*pstr == "Backward") {
        ident = BACKWARD_RECURRENCE_TIME;
      }
      else if (*pstr == "Forward") {
        ident = FORWARD_RECURRENCE_TIME;
      }
      else if (*pstr == "Mixture") {
        ident = MIXTURE;
      }
      else {
        genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                    "InterEvent or LengthBias or Backward or Forward or NbEvent or Mixture");
        return AMObj(AMObjType::ERROR);
      }
    }

    dist = ((Renewal*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , time);

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
    RWCString *pstr;
    bool status = true;
    int ident;
    Parametric_model *dist;
    Format_error error;


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
        ident = OBSERVATION;
      }
      else if (*pstr == "FirstOccurrence") {
        ident = FIRST_OCCURRENCE;
      }
      else if (*pstr == "Recurrence") {
        ident = RECURRENCE_TIME;
      }
      else if (*pstr == "Sojourn") {
        ident = SOJOURN_TIME;
      }
      else if (*pstr == "NbRun") {
        ident = NB_RUN;
      }
      else if (*pstr == "NbOccurrence") {
        ident = NB_OCCURRENCE;
      }
      else if (*pstr != "Forward") {
        status = false;
        genAMLError(ERRORMSG(DISTRIBUTION_NAME_sds) , "ExtractDistribution" , 2 ,
                    "Observation or FirstOccurrence or Recurrence or Sojourn or Forward or NbRun or NbOccurrence");
      }
    }

    if (*pstr == "Forward") {
      int nb_required = 3 , histogram_type = FINAL_RUN;


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
          if (*((AMString*)args[nb_required].val.p) != "HistogramType") {
            status = false;
            genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "ExtractDistribution" , nb_required + 1 ,
                        "HistogramType");
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
            histogram_type = INITIAL_RUN;
          }
          else if (*pstr == "FinalRun") {
            histogram_type = FINAL_RUN;
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
        dist = ((Semi_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[2].val.i , histogram_type);
        break;
      case AMObjType::HIDDEN_SEMI_MARKOV :
        dist = ((Hidden_semi_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[2].val.i , histogram_type);
        break;
      }
    }

    else {
      int variable , value;


      CHECKCONDVA((args.length() == 3) || (args.length() == 4) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractDistribution"));

      if ((ident == OBSERVATION) && (args.length() != 4)) {
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
        dist = ((Variable_order_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
        break;
      case AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV :
        dist = ((Hidden_variable_order_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
        break;
      case AMObjType::SEMI_MARKOV :
        dist = ((Semi_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
        break;
      case AMObjType::HIDDEN_SEMI_MARKOV :
        dist = ((Hidden_semi_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV :
        dist = ((Nonhomogeneous_markov*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , value);
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

  if (args[0].tag() == AMObjType::TOP_PARAMETERS) {
    Parametric_model *dist;
    Format_error error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractDistribution" , 2));

    CHECKCONDVA(args[1].tag() == AMObjType::INTEGER ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractDistribution" , 2 ,
                            args[1].tag.string().data() , "INT"));

    dist = ((Top_parameters*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[1].val.i);

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
 *  Extraction d'un histogramme elementaire,
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractHistogram(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

  if (args[0].tag() == AMObjType::MIXTURE_DATA) {
    RWCString *pstr;
    Distribution_data *histo;
    const Mixture_data *mixture_data;
    Format_error error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    mixture_data = (Mixture_data*)((STAT_model*)args[0].val.p)->pt;

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Weight") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      histo = new Distribution_data(*(mixture_data->get_weight()) ,
                                    mixture_data->get_mixture()->get_weight());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    if (*pstr == "Component") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
                              args[2].tag.string().data() , "INT"));

      histo = mixture_data->extract(error , args[2].val.i);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::HISTOGRAM , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Mixture") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      histo = new Distribution_data(*((Histogram*)mixture_data) , (Distribution*)mixture_data->get_mixture());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                "Weight or Component or Mixture");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::CONVOLUTION_DATA) {
    RWCString *pstr;
    Distribution_data *histo;
    const Convolution_data *convolution_data;
    Format_error error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    convolution_data = (Convolution_data*)((STAT_model*)args[0].val.p)->pt;

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Elementary") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
                            args[2].tag.string().data() , "INT"));

      histo = convolution_data->extract(error , args[2].val.i);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::HISTOGRAM , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
        return AMObj(AMObjType::ERROR);
      }
    }

    if (*pstr == "Convolution") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      histo = new Distribution_data(*((Histogram*)convolution_data) , (Distribution*)convolution_data->get_convolution());

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                "Elementary or Convolution");
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::COMPOUND_DATA) {
    RWCString *pstr;
    char type;
    Distribution_data *histo;
    Format_error error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Sum") {
      type = 's';
    }
    else if (*pstr == "Elementary") {
      type = 'e';
    }
    else if (*pstr == "Compound") {
      type = 'c';
    }
    else {
      genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                  "Sum or Elementary or Compound");
      return AMObj(AMObjType::ERROR);
    }

    histo = ((Compound_data*)((STAT_model*)args[0].val.p)->pt)->extract(error , type);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::VECTORS) {
    int nb_variable , variable;
    Distribution_data *histo;
    const Vectors *vec;
    Format_error error;


    vec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = vec->get_nb_variable();

    if (nb_variable == 1) {
      CHECKCONDVA(args.length() == 1 ,
                  genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ExtractHistogram"));
      variable = 1;
    }

    else {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      if (args[1].tag() != AMObjType::INTEGER) {
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
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
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::TIME_EVENTS) || (args[0].tag() == AMObjType::RENEWAL_DATA)) {
    RWCString *pstr;
    int ident;
    Distribution_data *histo;
    const Time_events *time_events;
    const Renewal_data *renewal_data;
    Format_error error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    time_events = (Time_events*)((STAT_model*)args[0].val.p)->pt;

    if (args[0].tag() == AMObjType::RENEWAL_DATA) {
      renewal_data = (Renewal_data*)((STAT_model*)args[0].val.p)->pt;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "NbEvent") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
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
        return AMObj(AMObjType::HISTOGRAM , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
        return AMObj(AMObjType::ERROR);
      }
    }

    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

    if (*pstr == "Mixture") {
      switch (args[0].tag()) {
      case AMObjType::TIME_EVENTS :
        histo = time_events->extract(error , MIXTURE);
        break;
      case AMObjType::RENEWAL_DATA :
        histo = renewal_data->extract(error , MIXTURE);
        break;
      }

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    if (*pstr == "ObservationTime") {
      histo = new Distribution_data(*(time_events->get_htime()));

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    if (args[0].tag() == AMObjType::TIME_EVENTS) {
      genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                  "ObservationTime or NbEvent or Mixture");
      return AMObj(AMObjType::ERROR);
    }

    CHECKCONDVA(args[0].tag() == AMObjType::RENEWAL_DATA ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 1 ,
                            args[0].tag.string().data() , "RENEWAL_DATA"));

    if (*pstr == "InterEvent") {
      ident = INTER_EVENT;
    }
    else if (*pstr == "Within") {
      ident = WITHIN_OBSERVATION_PERIOD;
    }
    else if (*pstr == "LengthBias") {
      ident = LENGTH_BIAS;
    }
    else if (*pstr == "Backward") {
      ident = BACKWARD_RECURRENCE_TIME;
    }
    else if (*pstr == "Forward") {
      ident = FORWARD_RECURRENCE_TIME;
    }
    else {
      genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                  "ObservationTime or InterEvent or Within or LengthBias or Backward or Forward or NbEvent or Mixture");
      return AMObj(AMObjType::ERROR);
    }

    histo = renewal_data->extract(error , ident);

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::SEQUENCES) || (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    RWCString *pstr;
    bool status = true;
    int ident , nb_variable , variable , offset , value;
    Distribution_data *histo;
    const Sequences *seq;
    const Markovian_sequences *markovian_seq;
    Format_error error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      seq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      seq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      seq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      seq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      seq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Length") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      histo = new Distribution_data(*(seq->get_hlength()));

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    if (*pstr == "Value") {
      nb_variable = seq->get_nb_variable();

      if (nb_variable == 1) {
        CHECKCONDVA(args.length() == 2 ,
                    genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));
        variable = 1;
      }

      else {
        CHECKCONDVA(args.length() == 3 ,
                    genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 3));

        if (args[2].tag() != AMObjType::INTEGER) {
          genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
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
        return AMObj(AMObjType::HISTOGRAM , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
        return AMObj(AMObjType::ERROR);
      }
    }

    CHECKCONDVA(args.length() >= 3 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      markovian_seq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      markovian_seq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      markovian_seq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      markovian_seq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    default :
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 1 , args[0].tag.string().data() ,
                  "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
      return AMObj(AMObjType::ERROR);
    }

    if (*pstr == "Observation") {
      ident = OBSERVATION;
    }
    else if (*pstr == "FirstOccurrence") {
      ident = FIRST_OCCURRENCE;
    }
    else if (*pstr == "Recurrence") {
      ident = RECURRENCE_TIME;
    }
    else if (*pstr == "Sojourn") {
      ident = SOJOURN_TIME;
    }
    else if (*pstr == "InitialRun") {
      ident = INITIAL_RUN;
    }
    else if (*pstr == "FinalRun") {
      ident = FINAL_RUN;
    }
    else if (*pstr == "NbRun") {
      ident = NB_RUN;
    }
    else if (*pstr == "NbOccurrence") {
      ident = NB_OCCURRENCE;
    }
    else {
      status = false;
      genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 ,
                  "Length or Observation or FirstOccurrence or Recurrence or Sojourn or InitialRun or FinalRun or NbRun or NbOccurrence");
    }

    if ((ident == OBSERVATION) && (args[0].tag() != AMObjType::VARIABLE_ORDER_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::SEMI_MARKOV_DATA) &&
        (args[0].tag() != AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 1 ,
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
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
                    args[2].tag.string().data() , "INT");
      }
      else {
        variable = args[2].val.i;
      }
    }

    CHECKCONDVA(args.length() == offset + 1 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , offset + 1));

    if (args[offset].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , offset + 1 ,
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
      histo = ((Markovian_sequences*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      histo = ((Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      histo = ((Semi_markov_data*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , variable , value);
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      histo = ((Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt)->extract(error , ident , value);
      break;
    }

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::TOPS) {
    RWCString *pstr;
    Distribution_data *histo;
    Format_error error;


    CHECKCONDVA(args.length() >= 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractHistogram"));

    CHECKCONDVA(args[1].tag() == AMObjType::STRING ,
                genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 2 ,
                            args[1].tag.string().data() , "STRING"));

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Main") {
      CHECKCONDVA(args.length() == 2 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 2));

      histo = new Distribution_data(*(((Tops*)((STAT_model*)args[0].val.p)->pt)->get_hlength()));

      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }

    if (*pstr == "Axillary") {
      CHECKCONDVA(args.length() == 3 ,
                  genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "ExtractHistogram" , 3));

      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "ExtractHistogram" , 3 ,
                              args[2].tag.string().data() , "INT"));

      histo = ((Tops*)((STAT_model*)args[0].val.p)->pt)->extract(error , args[2].val.i);

      if (histo) {
        STAT_model* model = new STAT_model(histo);
        return AMObj(AMObjType::HISTOGRAM , model);
      }
      else {
        AMLOUTPUT << "\n" << error;
        genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractHistogram");
        return AMObj(AMObjType::ERROR);
      }
    }

    genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractHistogram" , 2 , "Main or Axillary");
    return AMObj(AMObjType::ERROR);
  }

  genAMLError(ERRORMSG(STAT_DATA_TYPE_HISTOGRAM_sds) , "ExtractHistogram" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Extraction de mesures globales (longueur, temps avant la 1ere occurrence
 *  d'une valeur, nombre de series/d'occurrences d'une valeur) par sequence,
 *
 *--------------------------------------------------------------*/

AMObj STAT_ExtractVectors(const AMObjVector &args)

{
  RWCString *pstr;
  bool status = true;
  int type , nb_variable , variable , offset , value;
  Vectors *vec;
  const Sequences *seq;
  Format_error error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ExtractVectors"));

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    seq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    seq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    seq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    seq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    seq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

    vec = seq->extract_vectors(error , LENGTH);

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
    type = SEQUENCE_CUMUL;
  }
  else if (*pstr == "Mean") {
    type = SEQUENCE_MEAN;
  }
  else if (*pstr == "FirstOccurrence") {
    type = FIRST_OCCURRENCE;
  }
  else if (*pstr == "NbRun") {
    type = NB_RUN;
  }
  else if (*pstr == "NbOccurrence") {
    type = NB_OCCURRENCE;
  }
  else {
    status = false;
    genAMLError(ERRORMSG(HISTOGRAM_NAME_sds) , "ExtractVectors" , 2 ,
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

    vec = seq->extract_vectors(error , type , variable);

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

  vec = seq->extract_vectors(error , type , variable , value);

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

  case AMObjType::MIXTURE : {
    Mixture_data *mixt_histo;
    Format_error error;


    mixt_histo = ((Mixture*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

    if (mixt_histo) {
      STAT_model* model = new STAT_model(mixt_histo);
      return AMObj(AMObjType::MIXTURE_DATA , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ExtractData");
      return AMObj(AMObjType::ERROR);
    }
  }

  case AMObjType::CONVOLUTION : {
    Convolution_data *convol_histo;
    Format_error error;


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
    Compound_data *compound_histo;
    Format_error error;


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

  case AMObjType::HIDDEN_VARIABLE_ORDER_MARKOV : {
    Variable_order_markov_data *seq;
    Format_error error;


    seq = ((Hidden_variable_order_markov*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

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
    Semi_markov_data *seq;
    Format_error error;


    seq = ((Hidden_semi_markov*)((STAT_model*)args[0].val.p)->pt)->extract_data(error);

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
                "MIXTURE or CONVOLUTION or COMPOUND or HIDDEN_VARIABLE_ORDER_MARKOV or HIDDEN_SEMI-MARKOV");
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
  register int i;
  int nb_sample = args.length();


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Merge"));

  if ((args[0].tag() == AMObjType::HISTOGRAM) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    const Histogram **phisto;
    Distribution_data *histo;


    phisto = new const Histogram*[nb_sample];

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::HISTOGRAM :
        phisto[i] = (Histogram*)((Distribution_data*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::MIXTURE_DATA :
        phisto[i] = (Histogram*)((Mixture_data*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::CONVOLUTION_DATA :
        phisto[i] = (Histogram*)((Convolution_data*)((STAT_model*)args[i].val.p)->pt);
        break;
      case AMObjType::COMPOUND_DATA :
        phisto[i] = (Histogram*)((Compound_data*)((STAT_model*)args[i].val.p)->pt);
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 , args[i].tag.string().data() ,
                    "HISTOGRAM or MIXTURE_DATA or CONVOLUTION_DATA or COMPOUND_DATA");
        break;
      }
    }

    if (status) {
      histo = new Distribution_data(nb_sample , phisto);
    }

    delete [] phisto;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(histo);
    return AMObj(AMObjType::HISTOGRAM , model);
  }

  if (args[0].tag() == AMObjType::VECTORS) {
    const Vectors **pvec;
    Vectors *vec;
    Format_error error;


    pvec = new const Vectors*[nb_sample];
    pvec[0] = (Vectors*)((STAT_model*)args[0].val.p)->pt;

    for (i = 1;i < nb_sample;i++) {
      if (args[i].tag() != AMObjType::VECTORS) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "VECTORS");
      }
      else {
        pvec[i] = (Vectors*)((STAT_model*)args[i].val.p)->pt;
      }
    }

    if (status) {
      vec = pvec[0]->merge(error , nb_sample - 1 , pvec + 1);
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
    const Time_events **ptimev;
    Time_events *timev;


    ptimev = new const Time_events*[nb_sample];

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::TIME_EVENTS :
        ptimev[i] = (Time_events*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::RENEWAL_DATA :
        ptimev[i] = (Renewal_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "TIME_EVENTS or RENEWAL_DATA");
        break;
      }
    }

    if (status) {
      timev = new Time_events(nb_sample , ptimev);
    }

    delete [] ptimev;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(timev);
    return AMObj(AMObjType::TIME_EVENTS , model);
  }

  if (args[0].tag() == AMObjType::RENEWAL_DATA) {
    const Renewal_data **ptimev;
    Renewal_data *timev;
    Format_error error;


    ptimev = new const Renewal_data*[nb_sample];

    for (i = 0;i < nb_sample;i++) {
      if (args[i].tag() != AMObjType::RENEWAL_DATA) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "RENEWAL_DATA");
      }
      else {
        ptimev[i] = (Renewal_data*)((STAT_model*)args[i].val.p)->pt;
      }
    }

    if (status) {
      timev = ptimev[0]->merge(error , nb_sample - 1 , ptimev + 1);
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
    Sequences *seq;
    const Markovian_sequences **pmarkovian_seq;
    Markovian_sequences *markovian_seq;
    Format_error error;


    pseq = new const Sequences*[nb_sample];

    pmarkovian_seq = new const Markovian_sequences*[nb_sample];
    for (i = 0;i < nb_sample;i++) {
      pmarkovian_seq[i] = 0;
    }

    for (i = 0;i < nb_sample;i++) {
      switch (args[i].tag()) {
      case AMObjType::SEQUENCES :
        pseq[i] = (Sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MARKOVIAN_SEQUENCES :
        pseq[i] = (Markovian_sequences*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Markovian_sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
        pseq[i] = (Variable_order_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Variable_order_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::SEMI_MARKOV_DATA :
        pseq[i] = (Semi_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Semi_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
        pseq[i] = (Nonhomogeneous_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Nonhomogeneous_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 , args[i].tag.string().data() ,
                    "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
        break;
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
      seq = pseq[0]->merge(error , nb_sample - 1 , pseq + 1);
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
      markovian_seq = pmarkovian_seq[0]->merge(error , nb_sample - 1 , pmarkovian_seq + 1);
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

  if (args[0].tag() == AMObjType::TOPS) {
    const Tops **ptops;
    Tops *tops;


    ptops = new const Tops*[nb_sample];
    ptops[0] = (Tops*)((STAT_model*)args[0].val.p)->pt;

    for (i = 1;i < nb_sample;i++) {
      if (args[i].tag() != AMObjType::TOPS) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Merge" , i + 1 ,
                    args[i].tag.string().data() , "TOPS");
      }
      else {
        ptops[i] = (Tops*)((STAT_model*)args[i].val.p)->pt;
      }
    }

    if (status) {
      tops = new Tops(nb_sample , ptops);
    }

    delete [] ptops;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    STAT_model* model = new STAT_model(tops);
    return AMObj(AMObjType::TOPS , model);
  }

  if (args[0].tag() == AMObjType::CORRELATION) {
    const Correlation **pcorrel;
    Correlation *correl;
    Format_error error;


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

  if ((args[0].tag() == AMObjType::HISTOGRAM) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    bool status = true;
    const Histogram *ihisto;
    Distribution_data *histo;
    Format_error error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Shift" , 2));

    switch (args[0].tag()) {
    case AMObjType::HISTOGRAM :
      ihisto = (Histogram*)((Distribution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::MIXTURE_DATA :
      ihisto = (Histogram*)((Mixture_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (Histogram*)((Convolution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (Histogram*)((Compound_data*)((STAT_model*)args[0].val.p)->pt);
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
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Shift");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::SEQUENCES) ||
      (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    bool status = true;
    int nb_variable , variable , offset;
    const Vectors *ivec;
    const Sequences *iseq;
    Format_error error;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

    if (args[0].tag() == AMObjType::VECTORS) {
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
      Markovian_sequences *markovian_seq;


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

  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Shift" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
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

  if ((args[0].tag() == AMObjType::HISTOGRAM) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    RWCString *pstr;
    bool status = true;
    const Histogram *ihisto;
    Distribution_data *histo;
    Format_error error;


    CHECKCONDVA(args.length() == 3 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , 3));

    switch (args[0].tag()) {
    case AMObjType::HISTOGRAM :
      ihisto = (Histogram*)((Distribution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::MIXTURE_DATA :
      ihisto = (Histogram*)((Mixture_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (Histogram*)((Convolution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (Histogram*)((Compound_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    pstr = (AMString*)args[1].val.p;

    if (*pstr == "Step") {
      CHECKCONDVA(args[2].tag() == AMObjType::INTEGER ,
                  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Cluster" , 3 ,
                            args[2].tag.string().data() , "INT"));

      histo = ihisto->cluster(error , args[2].val.i);
    }

    else if (*pstr == "Information") {
      double ratio;


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

      histo = ihisto->cluster(error , ratio , AMLOUTPUT);
    }

    else if (*pstr == "Limit") {
      int nb_class = I_DEFAULT , *limit = 0;


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
      }

      if (status) {
        histo = ihisto->cluster(error , nb_class + 1 , limit);
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
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Cluster");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::SEQUENCES) ||
      (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    RWCString *pstr;
    bool status = true;
    int nb_required , nb_variable , variable , offset , nb_class = I_DEFAULT , *int_limit = 0;
    double *real_limit = 0;
    const Vectors *ivec;
    const Sequences *iseq;
    const Markovian_sequences *imarkovian_seq;
    Format_error error;


    // arguments obligatoires

    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      imarkovian_seq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      imarkovian_seq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      imarkovian_seq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      imarkovian_seq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
        }
        else if (parray->surfaceType() == AMObjType::REAL) {
          real_limit = buildRealArray(args , offset , "Cluster" , offset + 1 , nb_class , false);
          if (!real_limit) {
            status = false;
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

    if (args[0].tag() == AMObjType::VECTORS) {
      Vectors *vec;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
      }

      if (status) {
        if (*pstr == "Step") {
          vec = ivec->cluster(error , variable , args[offset].val.i);
        }

        else if (*pstr == "Limit") {
          if (int_limit) {
            vec = ivec->cluster(error , variable , nb_class + 1 , int_limit);
          }
          else if (real_limit) {
            vec = ivec->cluster(error , variable , nb_class + 1 , real_limit);
          }
        }
      }

      delete [] int_limit;
      delete [] real_limit;

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
      Markovian_sequences *markovian_seq;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
      }

      if (status) {
        if (*pstr == "Step") {
          seq = iseq->cluster(error , variable , args[offset].val.i);
        }

        else if (*pstr == "Limit") {
          if (int_limit) {
            seq = iseq->cluster(error , variable , nb_class + 1 , int_limit);
          }
          else if (real_limit) {
            seq = iseq->cluster(error , variable , nb_class + 1 , real_limit);
          }
        }
      }

      delete [] int_limit;
      delete [] real_limit;

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
      bool add_flag = false;
      Markovian_sequences *markovian_seq;


      pstr = (AMString*)args[1].val.p;

      if (*pstr == "Step") {
        if (args.length() != nb_required) {
          status = false;
          genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Cluster" , nb_required);
        }

        if (status) {
          markovian_seq = imarkovian_seq->cluster(error , variable , args[offset].val.i);
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
            add_flag = args[nb_required + 1].val.b;
          }

          if (real_limit) {
            status = false;
            genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "Cluster" , "AddVariable");
          }
        }

        if (status) {
          if (int_limit) {
            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
                                                    int_limit , add_flag);
          }
          else if (real_limit) {
            markovian_seq = imarkovian_seq->cluster(error , variable , nb_class + 1 ,
                                                    real_limit);
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

  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Cluster" , 1 ,
              args[0].tag.string().data());
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Transcodage des symboles.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Transcode(const AMObjVector &args)

{
  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Transcode"));

  if ((args[0].tag() == AMObjType::HISTOGRAM) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    bool status = true;
    int nb_symbol , *symbol = 0;
    const Histogram *ihisto;
    Distribution_data *histo;
    Format_error error;


    CHECKCONDVA(args.length() == 2 ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , 2));

    switch (args[0].tag()) {
    case AMObjType::HISTOGRAM :
      ihisto = (Histogram*)((Distribution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::MIXTURE_DATA :
      ihisto = (Histogram*)((Mixture_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (Histogram*)((Convolution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (Histogram*)((Compound_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    }

    nb_symbol = ihisto->nb_value - ihisto->offset;

    if (args[1].tag() != AMObjType::ARRAY) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , 2 ,
                  args[1].tag.string().data() , "ARRAY");
    }
    else {
      symbol = buildIntArray(args , 1 , "Transcode" , 2 , nb_symbol);
      if (!symbol) {
        status = false;
      }
    }

    if (status) {
      histo = ihisto->transcode(error , symbol);
    }

    delete [] symbol;

    if (!status) {
      return AMObj(AMObjType::ERROR);
    }

    if (histo) {
      STAT_model* model = new STAT_model(histo);
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Transcode");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::SEQUENCES) ||
      (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
      (args[0].tag() == AMObjType::VARIABLE_ORDER_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::SEMI_MARKOV_DATA) ||
      (args[0].tag() == AMObjType::NONHOMOGENEOUS_MARKOV_DATA)) {
    bool status = true;
    int nb_required , nb_variable , variable = I_DEFAULT , offset ,
        nb_symbol = I_DEFAULT , *symbol = 0;
    const Vectors *ivec;
    const Sequences *iseq;
    const Markovian_sequences *imarkovian_seq;
    Format_error error;


    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = ivec->get_nb_variable();
      break;
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      imarkovian_seq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      imarkovian_seq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      imarkovian_seq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = imarkovian_seq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      imarkovian_seq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

    switch (args[0].tag()) {
    case AMObjType::VECTORS :
      nb_symbol = (int)(ivec->get_max_value(variable - 1) - ivec->get_min_value(variable - 1)) + 1;
      break;
    case AMObjType::SEQUENCES :
      nb_symbol = (int)(iseq->get_max_value(variable - 1) - iseq->get_min_value(variable - 1)) + 1;
      break;
    default :
      nb_symbol = (int)(imarkovian_seq->get_max_value(variable - 1) - imarkovian_seq->get_min_value(variable - 1)) + 1;
      break;
    }

    if (args[offset].tag() != AMObjType::ARRAY) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Transcode" , offset + 1 ,
                    args[offset].tag.string().data() , "ARRAY");
    }
    else {
      symbol = buildIntArray(args , offset , "Transcode" , offset + 1 , nb_symbol);
      if (!symbol) {
        status = false;
      }
    }

    if (args[0].tag() == AMObjType::VECTORS) {
      Vectors *vec;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , nb_required);
      }

      if (status) {
        vec = ivec->transcode(error , variable , symbol);
      }

      delete [] symbol;

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
      Markovian_sequences *markovian_seq;


      if (args.length() != nb_required) {
        status = false;
        genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Transcode" , nb_required);
      }

      if (status) {
        seq = iseq->transcode(error , variable , symbol);
      }

      delete [] symbol;

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
      bool add_flag = false;
      Markovian_sequences *markovian_seq;


      if ((args.length() != nb_required) && (args.length() != nb_required + 2)) {
        delete [] symbol;
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
          add_flag = args[nb_required + 1].val.b;
        }
      }

      if (status) {
        markovian_seq = imarkovian_seq->transcode(error , variable , symbol , add_flag);
      }

      delete [] symbol;

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

  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "Transcode" , 1 ,
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
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , int_min_value , int_max_value;
  Format_error error;


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

  if ((args[0].tag() == AMObjType::HISTOGRAM) || (args[0].tag() == AMObjType::MIXTURE_DATA) ||
      (args[0].tag() == AMObjType::CONVOLUTION_DATA) || (args[0].tag() == AMObjType::COMPOUND_DATA)) {
    const Histogram *ihisto;
    Distribution_data *histo;


    CHECKCONDVA((nb_required == 2) || (nb_required == 3) ,
                genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "ValueSelect"));

    switch (args[0].tag()) {
    case AMObjType::HISTOGRAM :
      ihisto = (Histogram*)((Distribution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::MIXTURE_DATA :
      ihisto = (Histogram*)((Mixture_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::CONVOLUTION_DATA :
      ihisto = (Histogram*)((Convolution_data*)((STAT_model*)args[0].val.p)->pt);
      break;
    case AMObjType::COMPOUND_DATA :
      ihisto = (Histogram*)((Compound_data*)((STAT_model*)args[0].val.p)->pt);
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
      return AMObj(AMObjType::HISTOGRAM , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "ValueSelect");
      return AMObj(AMObjType::ERROR);
    }
  }

  if ((args[0].tag() == AMObjType::VECTORS) || (args[0].tag() == AMObjType::SEQUENCES) ||
      (args[0].tag() == AMObjType::MARKOVIAN_SEQUENCES) ||
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
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      nb_variable = iseq->get_nb_variable();
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

    if (args[0].tag() == AMObjType::VECTORS) {
      Vectors *vec;


      if (args[offset].tag() == AMObjType::INTEGER) {
        vec = ivec->value_select(error , variable , int_min_value , int_max_value , keep);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        vec = ivec->value_select(error , variable , real_min_value , real_max_value , keep);
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
      Markovian_sequences *markovian_seq;


      if (args[offset].tag() == AMObjType::INTEGER) {
        seq = iseq->value_select(error , variable , int_min_value , int_max_value , keep);
      }
      else if (args[offset].tag() == AMObjType::REAL) {
        seq = iseq->value_select(error , variable , real_min_value , real_max_value , keep);
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

  genAMLError(ERRORMSG(STAT_DATA_TYPE_TIME_EVENTS_TOPS_sds) , "ValueSelect" , 1 ,
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
  Format_error error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "VariableScaling"));

  switch (args[0].tag()) {
  case AMObjType::VECTORS :
    ivec = (Vectors*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = ivec->get_nb_variable();
    break;
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
    nb_variable = iseq->get_nb_variable();
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , 1 , args[0].tag.string().data() ,
                "VECTORS or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
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

  if (args[offset].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "VariableScaling" , offset + 1 ,
                args[offset].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::VECTORS) {
    Vectors *vec;


    vec = ivec->scaling(error , variable , args[offset].val.i);

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
    Markovian_sequences *markovian_seq;


    seq = iseq->scaling(error , variable , args[offset].val.i);

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
  bool status = true;
  int nb_required , variable = I_DEFAULT;
  Format_error error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Round"));

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Round" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Variable") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "Round" , nb_required + 1 , "Variable");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Round" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      variable = args[nb_required + 1].val.i;
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  if (args[0].tag() == AMObjType::VECTORS) {
    Vectors *vec;


    vec = ((Vectors*)((STAT_model*)args[0].val.p)->pt)->round(error , variable);

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
    Markovian_sequences *markovian_seq;


    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = iseq->round(error , variable);

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
              "VECTORS or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Selection par l'identificateur.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SelectIndividual(const AMObjVector &args)

{
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , nb_pattern = I_DEFAULT , *identifier = 0;
  Format_error error;


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

  if (args[0].tag() == AMObjType::VECTORS) {
    Vectors *vec;


    vec = ((Vectors*)((STAT_model*)args[0].val.p)->pt)->select_individual(error , nb_pattern ,
                                                                          identifier , keep);
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
    Markovian_sequences *markovian_seq;


    switch (args[0].tag()) {
    case AMObjType::SEQUENCES :
      iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = iseq->select_individual(error , nb_pattern , identifier , keep);
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

  if (args[0].tag() == AMObjType::TOPS) {
    Tops *tops;
    Format_error error;


    tops = ((Tops*)((STAT_model*)args[0].val.p)->pt)->select_individual(error , nb_pattern ,
                                                                        identifier , keep);
    delete [] identifier;

    if (tops) {
      STAT_model* model = new STAT_model(tops);
      return AMObj(AMObjType::TOPS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "SelectIndividual");
      return AMObj(AMObjType::ERROR);
    }
  }

  if (args[0].tag() == AMObjType::DISTANCE_MATRIX) {
    Distance_matrix *dist_matrix;


    dist_matrix = ((Distance_matrix*)((STAT_model*)args[0].val.p)->pt)->select_individual(error , nb_pattern ,
                                                                                          identifier , keep);
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
  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SelectIndividual" , 1 , args[0].tag.string().data() ,
              "VECTORS or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA or TOPS or DISTANCE_MATRIX");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Selection de variables.
 *
 *--------------------------------------------------------------*/

AMObj STAT_SelectVariable(const AMObjVector &args)

{
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , nb_variable = I_DEFAULT , *variable = 0;
  Format_error error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SelectVariable"));

  // argument obligatoire

  switch (args[1].tag()) {

  case AMObjType::INTEGER : {
    nb_variable = 1;
    variable = new int[nb_variable];
    variable[0] = args[1].val.i;
    break;
  }

  case AMObjType::ARRAY : {
    variable = buildIntArray(args , 1 , "SelectVariable" , 2 , nb_variable);
    if (!variable) {
      status = false;
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

  if (args[0].tag() == AMObjType::VECTORS) {
    Vectors *vec;


    vec = ((Vectors*)((STAT_model*)args[0].val.p)->pt)->select_variable(error , nb_variable ,
                                                                        variable , keep);
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
    Markovian_sequences *markovian_seq;


    seq = ((Sequences*)((STAT_model*)args[0].val.p)->pt)->select_variable(error , nb_variable ,
                                                                          variable , keep);
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
    const Markovian_sequences *iseq;
    Markovian_sequences *seq;


    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = iseq->select_variable(error , nb_variable , variable , keep);
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
              "VECTORS or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Concatenation des variables.
 *
 *--------------------------------------------------------------*/

AMObj STAT_MergeVariable(const AMObjVector &args)

{
  register int i;
  bool status = true;
  int nb_required , ref_sample = I_DEFAULT;
  Format_error error;


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

  if (args[0].tag() == AMObjType::VECTORS) {
    const Vectors **pvec;
    Vectors *vec;
    Format_error error;


    pvec = new const Vectors*[nb_required];
    pvec[0] = (Vectors*)((STAT_model*)args[0].val.p)->pt;

    for (i = 1;i < nb_required;i++) {
      if (args[i].tag() != AMObjType::VECTORS) {
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , i + 1 ,
                    args[i].tag.string().data() , "VECTORS");
      }
      else {
        pvec[i] = (Vectors*)((STAT_model*)args[i].val.p)->pt;
      }
    }

    if (status) {
      vec = pvec[0]->merge_variable(error , nb_required - 1 , pvec + 1 , ref_sample);
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
    Sequences *seq;
    const Markovian_sequences **pmarkovian_seq;
    Markovian_sequences *markovian_seq;


    pseq = new const Sequences*[nb_required];

    pmarkovian_seq = new const Markovian_sequences*[nb_required];
    for (i = 0;i < nb_required;i++) {
      pmarkovian_seq[i] = 0;
    }

    for (i = 0;i < nb_required;i++) {
      switch (args[i].tag()) {
      case AMObjType::SEQUENCES :
        pseq[i] = (Sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::MARKOVIAN_SEQUENCES :
        pseq[i] = (Markovian_sequences*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Markovian_sequences*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
        pseq[i] = (Variable_order_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Variable_order_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::SEMI_MARKOV_DATA :
        pseq[i] = (Semi_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Semi_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
        pseq[i] = (Nonhomogeneous_markov_data*)((STAT_model*)args[i].val.p)->pt;
        pmarkovian_seq[i] = (Nonhomogeneous_markov_data*)((STAT_model*)args[i].val.p)->pt;
        break;
      default :
        status = false;
        genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , i + 1 , args[i].tag.string().data() ,
                    "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
        break;
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
      seq = pseq[0]->merge_variable(error , nb_required - 1 , pseq + 1 , ref_sample);
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
      markovian_seq = pmarkovian_seq[0]->merge_variable(error , nb_required - 1 ,
                                                        pmarkovian_seq + 1 , ref_sample);
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

  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MergeVariable" , 1 , args[0].tag.string().data() ,
              "VECTORS or SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
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
  const Time_events *itimev;
  Time_events *timev;
  Format_error error;


  CHECKCONDVA(args.length() == 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "NbEventSelect" , 3));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (Time_events*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (Renewal_data*)((STAT_model*)args[0].val.p)->pt;
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
  const Time_events *itimev;
  Time_events *timev;
  Format_error error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "TimeScaling" , 2));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (Time_events*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (Renewal_data*)((STAT_model*)args[0].val.p)->pt;
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
  const Time_events *itimev;
  Time_events *timev;
  Format_error error;


  CHECKCONDVA((args.length() == 2) || (args.length() == 3) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "TimeSelect"));

  switch (args[0].tag()) {
  case AMObjType::TIME_EVENTS :
    itimev = (Time_events*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::RENEWAL_DATA :
    itimev = (Renewal_data*)((STAT_model*)args[0].val.p)->pt;
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
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , min_length , max_length;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


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
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

  seq = iseq->length_select(error , min_length , max_length , keep);

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
 *  Suppression des premieres/dernieres series d'une valeur donne.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RemoveRun(const AMObjVector &args)

{
  RWCString *pstr;
  char position;
  bool status = true;
  int nb_required , nb_variable , variable , offset , value , max_run_length = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


  CHECKCONDVA(args.length() >= 3 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RemoveRun"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
      position = 'b';
    }
    else if (*pstr == "End") {
      position = 'e';
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
    Format_error error;


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
    Markovian_sequences *iseq , *seq;


    switch (args[0].tag()) {
    case AMObjType::MARKOVIAN_SEQUENCES :
      iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
      iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::SEMI_MARKOV_DATA :
      iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
      iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
      break;
    }

    seq = new Markovian_sequences(*iseq , 'c' , REVERSE);

    STAT_model* model = new STAT_model(seq);
    return AMObj(AMObjType::MARKOVIAN_SEQUENCES , model);
  }

  if (args[0].tag() == AMObjType::TOPS) {
    Tops *tops;
    Format_error error;


    tops = ((Tops*)((STAT_model*)args[0].val.p)->pt)->reverse(error);

    if (tops) {
      STAT_model* model = new STAT_model(tops);
      return AMObj(AMObjType::TOPS , model);
    }
    else {
      AMLOUTPUT << "\n" << error;
      genAMLError(ERRORMSG(STAT_MODULE_s) , "Reverse");
      return AMObj(AMObjType::ERROR);
    }
  }

  genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "Reverse" , args[0].tag.string().data() ,
              "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA or TOPS");
  return AMObj(AMObjType::ERROR);
}


/*--------------------------------------------------------------*
 *
 *  Suppression du parametre d'index.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RemoveIndexParameter(const AMObjVector &args)

{
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "RemoveIndexParameter"));

  // argument obligatoire

  switch (args[0].tag()) {

  case AMObjType::SEQUENCES : {
    Sequences *seq;
    Markovian_sequences *markovian_seq;


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
    Markovian_sequences *seq;


    seq = ((Markovian_sequences*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

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
    Variable_order_markov_data *seq;


    seq = ((Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

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
    Semi_markov_data *seq;


    seq = ((Semi_markov_data*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

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
    Nonhomogeneous_markov_data *seq;


    seq = ((Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt)->remove_index_parameter(error);

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
 *  Selection sur les valeurs prises par le parametre d'index.
 *
 *--------------------------------------------------------------*/

AMObj STAT_IndexParameterSelect(const AMObjVector &args)

{
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , min_index_parameter , max_index_parameter;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


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
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

  seq = iseq->index_parameter_select(error , min_index_parameter , max_index_parameter , keep);

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
 *  Extraction de sous-sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_IndexParameterExtract(const AMObjVector &args)

{
  bool status = true;
  int nb_required , max_index_parameter = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "IndexParameterExtract"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  RWCString *pstr;
  bool status = true , keep = true;
  int nb_required , nb_value = I_DEFAULT , *value = 0;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


  nb_required = 3;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SegmentationExtract"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) != "Mode") {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "SegmentationExtract" , nb_required + 1 , "Mode");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::STRING) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "SegmentationExtract" , nb_required + 1 ,
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
        genAMLError(ERRORMSG(ARG_sds) , "SegmentationExtract" , nb_required + 1 , "Keep or Reject");
      }
    }
  }

  if (status) {
    seq = iseq->segmentation_extract(error , args[1].val.i , nb_value , value , keep);
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
  Markovian_sequences *markovian_seq;
  Format_error error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "Cumulate"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  RWCString *pstr;
  char *file_name = 0;
  bool status = true , variable_option = false , first_element_option = false ,
       first_element = false , file_name_option = false;
  register int i;
  int nb_required , variable = I_DEFAULT;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


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
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
        switch (variable_option) {

        case false : {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Difference" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Difference" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "FirstElement") {
        switch (first_element_option) {

        case false : {
          first_element_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Difference" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            first_element = args[nb_required + i * 2 + 1].val.b;
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "Difference" , nb_required + i + 1);
          break;
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
 *  Filtrage de type moyenne mobile des sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_MovingAverage(const AMObjVector &args)

{
  RWCString *pstr;
  bool status = true , variable_option = false , begin_end_option = false , begin_end = false ,
       output_option = false;
  register int i;
  int nb_required , nb_point = I_DEFAULT , int_sum , variable = I_DEFAULT ,
      output = TREND , *int_filter;
  double sum , *filter;
  const Distribution *dist;
  const Sequences *iseq;
  Sequences *seq;
  Markovian_sequences *markovian_seq;
  Format_error error;


  nb_required = 2;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) || (args.length() == nb_required + 6) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "MovingAverage"));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "MovingAverage" , args[0].tag.string().data() ,
                "SEQUENCES or MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    break;
  }

  filter = 0;
  dist = 0;

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
  }

  else {
    switch (args[1].tag()) {
    case AMObjType::DISTRIBUTION :
      dist = new Distribution(*((Distribution*)((Parametric_model*)((STAT_model*)args[1].val.p)->pt)));
      break;
    case AMObjType::MIXTURE :
      dist = new Distribution(*((Distribution*)((Mixture*)((STAT_model*)args[1].val.p)->pt)));
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
                  "ARRAY or DISTRIBUTION or MIXTURE or CONVOLUTION or COMPOUND");
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
        switch (variable_option) {

        case false : {
          variable_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::INTEGER) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "INT");
          }
          else {
            variable = args[nb_required + i * 2 + 1].val.i;
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "BeginEnd") {
        switch (begin_end_option) {

        case false : {
          begin_end_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "MovingAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            begin_end = args[nb_required + i * 2 + 1].val.b;
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "Output") {
        switch (output_option) {

        case false : {
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
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "MovingAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "MovingAverage" , nb_required + i + 1 ,
                    "Variable or BeginEnd or Output");
      }
    }
  }

  if (status) {
    if (!dist) {
      seq = iseq->moving_average(error , nb_point , filter , variable , begin_end , output);
    }
    else {
      seq = iseq->moving_average(error , *dist , variable , begin_end , output);
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
  RWCString *pstr;
  char format = 'a' , *file_name = 0;
  bool status = true , standard_deviation_option = false , standard_deviation = false ,
       output_option = false , file_name_option = false , format_option = false;
  register int i;
  int nb_required , output = SEQUENCE;
  const Sequences *iseq;
  Sequences *seq;
  Format_error error;


  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ||
              (args.length() == nb_required + 4) || (args.length() == nb_required + 6) ||
              (args.length() == nb_required + 8) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "PointwiseAverage"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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

      if (*pstr == "StandardDeviation") {
        switch (standard_deviation_option) {

        case false : {
          standard_deviation_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::BOOL) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "BOOL");
          }
          else {
            standard_deviation = args[nb_required + i * 2 + 1].val.b;
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "Output") {
        switch (output_option) {

        case false : {
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
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "FileName") {
        switch (file_name_option) {

        case false : {
          file_name_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            file_name = (char*)((AMString*)args[nb_required + i * 2 + 1].val.p)->data();
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else if (*pstr == "Format") {
        switch (format_option) {

        case false : {
          format_option = true;

          if (args[nb_required + i * 2 + 1].tag() != AMObjType::STRING) {
            status = false;
            genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "PointwiseAverage" , nb_required + i + 1 ,
                        args[nb_required + i * 2 + 1].tag.string().data() , "STRING");
          }
          else {
            pstr = (AMString*)args[nb_required + i * 2 + 1].val.p;
            if (*pstr == "ASCII") {
              format = 'a';
            }
            else if (*pstr == "SpreadSheet") {
              format = 's';
            }
            else {
              status = false;
              genAMLError(ERRORMSG(K_FILE_FORMAT_ERR_sds) , "PointwiseAverage" ,
                          nb_required + i + 1 , "ASCII or SpreadSheet");
            }
          }
          break;
        }

        case true : {
          status = false;
          genAMLError(ERRORMSG(USED_OPTION_sd) , "PointwiseAverage" , nb_required + i + 1);
          break;
        }
        }
      }

      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "PointwiseAverage" , nb_required + i + 1 ,
                    "StandardDeviation or Output or FileName or Format");
      }
    }
  }

  if ((format_option) && (!file_name_option)) {
    status = false;
    genAMLError(ERRORMSG(FORBIDDEN_OPTION_ss) , "PointwiseAverage" , "Format");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->pointwise_average(error , standard_deviation , output , file_name , format);

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
  Markovian_sequences *markovian_seq;
  Format_error error;


  CHECKCONDVA(args.length() >= 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "RecurrenceTimeSequences"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  Markovian_sequences *markovian_seq;
  Format_error error;


  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "SojournTimeSequences"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  Markovian_sequences *markovian_seq;
  Format_error error;


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
  Markovian_sequences *markovian_seq;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Cross"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::SEQUENCES :
    iseq = (Sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  const Markovian_sequences *iseq;
  Markovian_sequences *seq;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ComputeInitialRun"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  bool status = true;
  int nb_required , sequence_length = I_DEFAULT , run_length = I_DEFAULT;
  const Markovian_sequences *iseq;
  Markovian_sequences *seq;
  Format_error error;


  CHECKCONDVA(args.length() >= 1 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "AddAbsorbingRun"));

  // argument obligatoire

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  default :
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "AddAbsorbingRun" , args[0].tag.string().data() ,
                "MARKOVIAN_SEQUENCES or VARIABLE_ORDER_MARKOV_DATA or SEMI-MARKOV_DATA or NONHOMOGENEOUS_MARKOV_DATA");
    return AMObj(AMObjType::ERROR);
  }

  nb_required = 1;

  CHECKCONDVA((args.length() == nb_required) || (args.length() == nb_required + 2) ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_s) , "AddAbsorbingRun"));

  // argument optionnel

  if (args.length() == nb_required + 2) {
    if (args[nb_required].tag() != AMObjType::OPTION) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + 1 ,
                  args[nb_required].tag.string().data() , "OPTION");
    }
    else {
      if (*((AMString*)args[nb_required].val.p) == "SequenceLength") {
        sequence_length = 0;
      }
      else if (*((AMString*)args[nb_required].val.p) == "RunLength") {
        run_length = 0;
      }
      else {
        status = false;
        genAMLError(ERRORMSG(K_OPTION_NAME_ERR_sds) , "AddAbsorbingRun" , nb_required + 1 ,
                    "SequenceLength or RunLength");
      }
    }

    if (args[nb_required + 1].tag() != AMObjType::INTEGER) {
      status = false;
      genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "AddAbsorbingRun" , nb_required + 1 ,
                  args[nb_required + 1].tag.string().data() , "INT");
    }
    else {
      if (sequence_length == 0) {
        sequence_length = args[nb_required + 1].val.i;
      }
      else if (run_length == 0){
        run_length = args[nb_required + 1].val.i;
      }
    }
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  seq = iseq->add_absorbing_run(error , sequence_length , run_length);

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
 *  Segmentation de sequences.
 *
 *--------------------------------------------------------------*/

AMObj STAT_Split(const AMObjVector &args)

{
  bool status = true;
  const Markovian_sequences *iseq;
  Markovian_sequences *seq;
  Format_error error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "Split" , 2));

  // arguments obligatoires

  switch (args[0].tag()) {
  case AMObjType::MARKOVIAN_SEQUENCES :
    iseq = (Markovian_sequences*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::VARIABLE_ORDER_MARKOV_DATA :
    iseq = (Variable_order_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::SEMI_MARKOV_DATA :
    iseq = (Semi_markov_data*)((STAT_model*)args[0].val.p)->pt;
    break;
  case AMObjType::NONHOMOGENEOUS_MARKOV_DATA :
    iseq = (Nonhomogeneous_markov_data*)((STAT_model*)args[0].val.p)->pt;
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
  Distance_matrix *dist_matrix;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "ToDistanceMatrix"));

  CHECKCONDVA(args[0].tag() == AMObjType::CLUSTERS ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sss) , "ToDistanceMatrix" ,
                          args[0].tag.string().data() , "CLUSTERS"));

  dist_matrix = new Distance_matrix(*((Clusters*)((STAT_model*)args[0].val.p)->pt));

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
  Distance_matrix *dist_matrix;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Symmetrize"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTANCE_MATRIX ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Symmetrize" , 1 ,
                          args[0].tag.string().data() , "DISTANCE_MATRIX"));

  dist_matrix = ((Distance_matrix*)((STAT_model*)args[0].val.p)->pt)->symmetrize(error);

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
  Distance_matrix *dist_matrix;
  Format_error error;


  CHECKCONDVA(args.length() == 1 ,
              genAMLError(ERRORMSG(K_SINGLE_ARG_ERR_s) , "Unnormalize"));

  CHECKCONDVA(args[0].tag() == AMObjType::DISTANCE_MATRIX ,
              genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "Unnormalize" , 1 ,
                          args[0].tag.string().data() , "DISTANCE_MATRIX"));

  dist_matrix = ((Distance_matrix*)((STAT_model*)args[0].val.p)->pt)->unnormalize(error);

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


/*--------------------------------------------------------------*
 *
 *  Suppression des premiers entrenoeuds de l'axe porteur.
 *
 *--------------------------------------------------------------*/

AMObj STAT_RemoveApicalInternodes(const AMObjVector &args)

{
  bool status = true;
  Tops *tops;
  Format_error error;


  CHECKCONDVA(args.length() == 2 ,
              genAMLError(ERRORMSG(K_NB_ARG_ERR_sd) , "RemoveApicalInternodes" , 2));

  if (args[0].tag() != AMObjType::TOPS) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveApicalInternodes" , 1 ,
                args[0].tag.string().data() , "TOPS");
  }

  if (args[1].tag() != AMObjType::INTEGER) {
    status = false;
    genAMLError(ERRORMSG(K_F_ARG_TYPE_ERR_sdss) , "RemoveApicalInternodes" , 2 ,
                args[1].tag.string().data() , "INT");
  }

  if (!status) {
    return AMObj(AMObjType::ERROR);
  }

  tops = ((Tops*)((STAT_model*)args[0].val.p)->pt)->shift(error , args[1].val.i);

  if (tops) {
    STAT_model* model = new STAT_model(tops);
    return AMObj(AMObjType::TOPS , model);
  }
  else {
    AMLOUTPUT << "\n" << error;
    genAMLError(ERRORMSG(STAT_MODULE_s) , "RemoveApicalInternodes");
    return AMObj(AMObjType::ERROR);
  }
}
