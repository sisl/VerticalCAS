/*
 ------------------------------------------------------------------
 ** Top contributors:
 **   Shiqi Wang and Suman Jana
 ** This file is part of the ReluVal project.
 ** Copyright (c) 2018-2019 by the authors listed in the file LICENSE
 ** and their institutional affiliations.
 ** All rights reserved.
 -----------------------------------------------------------------
 *
 * This is the main file of ReluVal, here is the usage:
 * ./network_test [property] [network] [target] 
 *      [need to print=0] [test for one run=0] [check mode=0]
 *
 * [property]: the saftety property want to verify
 *
 * [network]: the network want to test with
 *
 * [target]: Wanted label of the property
 *
 * [need to print]: whether need to print the detailed info of each split.
 * 0 is not and 1 is yes. Default value is 0.
 *
 * [test for one run]: whether need to estimate the output range without
 * split. 0 is no, 1 is yes. Default value is 0.
 *
 * [check mode]: normal split mode is 0. Check adv mode is 1.
 * Check adv mode will prevent further splits as long as the depth goes
 * upper than 20 so as to locate the concrete adversarial examples faster.
 * Default value is 0.
 * 
 */


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "split_VC.h"

//extern int thread_tot_cnt;

/* print the progress if getting SIGQUIT */
void sig_handler(int signo)
{

    if (signo == SIGQUIT) {
        printf("progress: %d/1024\n", progress);
    }

}



int main( int argc, char *argv[])
{

    //char *FULL_NET_PATH =\
            "nnet/ACASXU_run2a_1_1_batch_2000.nnet";
    char *FULL_NET_PATH;

    PROPERTY=2;

    if (argc > 9 || argc < 2) {
        printf("please specify a network\n");
        printf("./network_test [network] "
            "[print] "
            "[test for one run] [check mode]\n");
        exit(1);
    }

    for (int i=1;i<argc;i++) {

        if (i == 1) {
            FULL_NET_PATH = argv[i];
        }

        if (i == 2) {
            NEED_PRINT = atoi(argv[i]);
            if(NEED_PRINT != 0 && NEED_PRINT!=1){
                printf("Wrong print");
                exit(1);
            }
        }

        if (i == 3) {
            NEED_FOR_ONE_RUN = atoi(argv[i]);

            if (NEED_FOR_ONE_RUN != 0 && NEED_FOR_ONE_RUN != 1) {
                printf("Wrong test for one run");
                exit(1);
            }

        }

        if (i == 4) {

            if (atoi(argv[i]) == 0) {
                CHECK_ADV_MODE = 0;
                PARTIAL_MODE = 0;
            }

            if (atoi(argv[i]) == 1) {
                CHECK_ADV_MODE = 1;
                PARTIAL_MODE = 0;
            }

            if (atoi(argv[i]) == 2) {
                CHECK_ADV_MODE = 0;
                PARTIAL_MODE = 1;
            }

        }

    }

    openblas_set_num_threads(1);

    
    
    //clock_t start, end;
    srand((unsigned)time(NULL));
    double time_spent;
    int i,j,layer;

    
    struct NNet* nnet = load_network(FULL_NET_PATH, 0);
    int numLayers    = nnet->numLayers;
    int inputSize    = nnet->inputSize;
    int outputSize   = nnet->outputSize;
    
    float u[inputSize], l[inputSize];

    int n = 0;
    int feature_range_length = 0;
    int split_feature = -1;
    int depth = 0;

    float vints[]  = {-25.000,-24.000,-23.000,-22.000,-21.000,-20.000,-19.000,-18.000,-17.000,-16.000,-15.000,-14.000,-13.000,-12.000,-11.000,-10.000,-9.000,-8.000,-7.000,-6.000,-5.000,-4.000,-3.000,-2.000,-1.000,0.000,1.000,2.000,3.000,4.000,5.000,6.000,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,17.000,18.000,19.000,20.000,21.000,22.000,23.000,24.000,25.000};
    float vowns[]  = {-41.667,-41.666,-41.000,-40.000,-39.000,-38.000,-37.000,-36.000,-35.000,-34.000,-33.000,-32.000,-31.000,-30.000,-29.000,-28.000,-27.000,-26.000,-25.001,-24.999,-24.000,-23.000,-22.000,-21.000,-20.000,-19.000,-18.000,-17.000,-16.000,-15.000,-14.000,-13.000,-12.000,-11.000,-10.000,-9.000,-8.000,-7.000,-6.000,-5.000,-4.000,-3.000,-2.000,-1.000,-0.001,0.001,1.000,2.000,3.000,4.000,5.000,6.000,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,17.000,18.000,19.000,20.000,21.000,22.000,23.000,24.000,24.999,25.001,26.000,27.000,28.000,29.000,30.000,31.000,32.000,33.000,34.000,35.000,36.000,37.000,38.000,39.000,40.000,41.000,41.666,41.667};
    float hs[] = {-3000.000,-2900.000,-2800.000,-2700.000,-2600.000,-2500.000,-2400.000,-2300.000,-2200.000,-2100.000,-2000.000,-1900.000,-1800.000,-1700.000,-1600.000,-1500.000,-1450.000,-1400.000,-1350.000,-1300.000,-1250.000,-1200.000,-1150.000,-1100.000,-1050.000,-1000.000,-950.000,-900.000,-875.000,-850.000,-825.000,-800.000,-775.000,-750.000,-725.000,-700.000,-690.000,-680.000,-670.000,-660.000,-650.000,-640.000,-630.000,-620.000,-610.000,-600.000,-590.000,-580.000,-570.000,-560.000,-550.000,-540.000,-530.000,-520.000,-510.000,-500.000,-497.000,-494.000,-491.000,-488.000,-485.000,-482.000,-479.000,-476.000,-473.000,-470.000,-467.000,-464.000,-461.000,-458.000,-455.000,-452.000,-449.000,-446.000,-443.000,-440.000,-437.000,-434.000,-431.000,-428.000,-425.000,-422.000,-419.000,-416.000,-413.000,-410.000,-407.000,-404.000,-401.000,-398.000,-395.000,-392.000,-389.000,-386.000,-383.000,-380.000,-377.000,-374.000,-371.000,-368.000,-365.000,-362.000,-359.000,-356.000,-353.000,-350.000,-347.000,-344.000,-341.000,-338.000,-335.000,-332.000,-329.000,-326.000,-323.000,-320.000,-317.000,-314.000,-311.000,-308.000,-305.000,-302.000,-299.000,-296.000,-293.000,-290.000,-287.000,-284.000,-281.000,-278.000,-275.000,-272.000,-269.000,-266.000,-263.000,-260.000,-257.000,-254.000,-251.000,-248.000,-245.000,-242.000,-239.000,-236.000,-233.000,-230.000,-227.000,-224.000,-221.000,-218.000,-215.000,-212.000,-209.000,-206.000,-203.000,-200.000,-198.000,-196.000,-194.000,-192.000,-190.000,-188.000,-186.000,-184.000,-182.000,-180.000,-178.000,-176.000,-174.000,-172.000,-170.000,-168.000,-166.000,-164.000,-162.000,-160.000,-158.000,-156.000,-154.000,-152.000,-150.000,-148.000,-146.000,-144.000,-142.000,-140.000,-138.000,-136.000,-134.000,-132.000,-130.000,-128.000,-126.000,-124.000,-122.000,-120.000,-118.000,-116.000,-114.000,-112.000,-110.000,-108.000,-106.000,-104.000,-102.000,-100.000,-99.000,-98.000,-97.000,-96.000,-95.000,-94.000,-93.000,-92.000,-91.000,-90.000,-89.000,-88.000,-87.000,-86.000,-85.000,-84.000,-83.000,-82.000,-81.000,-80.000,-79.000,-78.000,-77.000,-76.000,-75.000,-74.000,-73.000,-72.000,-71.000,-70.000,-69.000,-68.000,-67.000,-66.000,-65.000,-64.000,-63.000,-62.000,-61.000,-60.000,-59.000,-58.000,-57.000,-56.000,-55.000,-54.000,-53.000,-52.000,-51.000,-50.000,-49.000,-48.000,-47.000,-46.000,-45.000,-44.000,-43.000,-42.000,-41.000,-40.000,-39.000,-38.000,-37.000,-36.000,-35.000,-34.000,-33.000,-32.000,-31.000,-30.000,-29.000,-28.000,-27.000,-26.000,-25.000,-24.000,-23.000,-22.000,-21.000,-20.000,-19.000,-18.000,-17.000,-16.000,-15.000,-14.000,-13.000,-12.000,-11.000,-10.000,-9.000,-8.000,-7.000,-6.000,-5.000,-4.000,-3.000,-2.000,-1.000,0.000,1.000,2.000,3.000,4.000,5.000,6.000,7.000,8.000,9.000,10.000,11.000,12.000,13.000,14.000,15.000,16.000,17.000,18.000,19.000,20.000,21.000,22.000,23.000,24.000,25.000,26.000,27.000,28.000,29.000,30.000,31.000,32.000,33.000,34.000,35.000,36.000,37.000,38.000,39.000,40.000,41.000,42.000,43.000,44.000,45.000,46.000,47.000,48.000,49.000,50.000,51.000,52.000,53.000,54.000,55.000,56.000,57.000,58.000,59.000,60.000,61.000,62.000,63.000,64.000,65.000,66.000,67.000,68.000,69.000,70.000,71.000,72.000,73.000,74.000,75.000,76.000,77.000,78.000,79.000,80.000,81.000,82.000,83.000,84.000,85.000,86.000,87.000,88.000,89.000,90.000,91.000,92.000,93.000,94.000,95.000,96.000,97.000,98.000,99.000,100.000,102.000,104.000,106.000,108.000,110.000,112.000,114.000,116.000,118.000,120.000,122.000,124.000,126.000,128.000,130.000,132.000,134.000,136.000,138.000,140.000,142.000,144.000,146.000,148.000,150.000,152.000,154.000,156.000,158.000,160.000,162.000,164.000,166.000,168.000,170.000,172.000,174.000,176.000,178.000,180.000,182.000,184.000,186.000,188.000,190.000,192.000,194.000,196.000,198.000,200.000,203.000,206.000,209.000,212.000,215.000,218.000,221.000,224.000,227.000,230.000,233.000,236.000,239.000,242.000,245.000,248.000,251.000,254.000,257.000,260.000,263.000,266.000,269.000,272.000,275.000,278.000,281.000,284.000,287.000,290.000,293.000,296.000,299.000,302.000,305.000,308.000,311.000,314.000,317.000,320.000,323.000,326.000,329.000,332.000,335.000,338.000,341.000,344.000,347.000,350.000,353.000,356.000,359.000,362.000,365.000,368.000,371.000,374.000,377.000,380.000,383.000,386.000,389.000,392.000,395.000,398.000,401.000,404.000,407.000,410.000,413.000,416.000,419.000,422.000,425.000,428.000,431.000,434.000,437.000,440.000,443.000,446.000,449.000,452.000,455.000,458.000,461.000,464.000,467.000,470.000,473.000,476.000,479.000,482.000,485.000,488.000,491.000,494.000,497.000,500.000,510.000,520.000,530.000,540.000,550.000,560.000,570.000,580.000,590.000,600.000,610.000,620.000,630.000,640.000,650.000,660.000,670.000,680.000,690.000,700.000,725.000,750.000,775.000,800.000,825.000,850.000,875.000,900.000,950.000,1000.000,1050.000,1100.000,1150.000,1200.000,1250.000,1300.000,1350.000,1400.000,1450.000,1500.000,1600.000,1700.000,1800.000,1900.000,2000.000,2100.000,2200.000,2300.000,2400.000,2500.000,2600.000,2700.000,2800.000,2900.000,3000.000};
    
    
    int lenVINTS = 50;
    int lenVOWNS = 89;
    int lenHS = 610;
    

    printf("running with network %s\n", FULL_NET_PATH);

    
    for (int target=0; target<9; target++){
        struct NNet* nnet = load_network(FULL_NET_PATH, target);
        for (float tau=0.0; tau<40.5; tau+=1){
            printf("RA: %d, TAU: %.1f\n",target,tau);
            fflush(stdout); 
            
            gettimeofday(&start, NULL);
            for (int reg=0; reg<lenVINTS*lenVOWNS*lenHS; reg++){ //9890
                
                reset_variables();
                n = 0;
                feature_range_length = 0;
                split_feature = -1;
                depth = 0;


                load_inputs(reg, inputSize, u, l,tau, hs, vowns, vints, lenHS, lenVOWNS, lenVINTS);

                struct Matrix input_upper = {u,1,nnet->inputSize};
                struct Matrix input_lower = {l,1,nnet->inputSize};

                struct Interval input_interval = {input_lower, input_upper};

                float grad_upper[inputSize], grad_lower[inputSize];
                struct Interval grad_interval = {
                            (struct Matrix){grad_upper, 1, inputSize},
                            (struct Matrix){grad_lower, 1, inputSize}
                        };

                normalize_input_interval(nnet, &input_interval);

                float o[nnet->outputSize];
                struct Matrix output = {o, outputSize, 1};

                float o_upper[nnet->outputSize], o_lower[nnet->outputSize];
                struct Interval output_interval = {
                            (struct Matrix){o_lower, outputSize, 1},
                            (struct Matrix){o_upper, outputSize, 1}
                        };


                for (int i=0;i<inputSize;i++) {

                    if (input_interval.upper_matrix.data[i] <\
                            input_interval.lower_matrix.data[i]) {
                        printf("wrong input!\n");
                        printf("Index: %d, Upper: %.5e, Lower: %.5e\n",i,input_interval.upper_matrix.data[i],input_interval.lower_matrix.data[i]);
                        exit(0);
                    }

                    if(input_interval.upper_matrix.data[i] !=\
                            input_interval.lower_matrix.data[i]){
                        n++;
                    }

                }



                feature_range_length = n;
                int *feature_range = (int*)malloc(n*sizeof(int));

                for (int i=0, n=0;i<nnet->inputSize;i++) {
                    if(input_interval.upper_matrix.data[i] !=\
                            input_interval.lower_matrix.data[i]){
                        feature_range[n] = i;
                        n++;
                    }
                }

                int isOverlap = 0;
                float avg[100] = {0};
                adv_found = 0;
                for (int i=0;i<1;i++) {
                    //forward_prop_interval_equation(nnet,\
                            &input_interval, &output_interval,\
                            &grad_interval);
                    isOverlap = direct_run_check(nnet,\
                            &input_interval, &output_interval,\
                            &grad_interval, depth, feature_range,\
                            feature_range_length, split_feature);
                }




                if (isOverlap == 0 && adv_found == 0 && depth_exceeded==0) {
                    //printf("UNSAT: %d\n",reg);
                } else if (adv_found!=0) {
                    printf("%d,\n",reg);
                } else if (depth_exceeded==1){
                    printf("%d,\n",reg);
                    //printf("DEPTH EXCEEDED: ");
                    //printf("%d,\n",reg);
                }
                free(feature_range);
            }

            gettimeofday(&finish, NULL);
            time_spent = ((float)(finish.tv_sec - start.tv_sec) *\
                        1000000 + (float)(finish.tv_usec - start.tv_usec)) /\
                        1000000;

            printf("time: %f \n\n\n", time_spent);
        }
    }

    destroy_network(nnet);

}
