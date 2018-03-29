//
//  p105_data.c
//  X_euler
//
//  Created by Jeannot on 22/03/2018.
//  Copyright © 2018 Jeannot. All rights reserved.
//

#include "p105_data.h"
static int32_t subset00[] = { 81,88,75,42,87,84,86,65,0 } ;
static int32_t subset01[] = { 157,150,164,119,79,159,161,139,158,0 } ;
static int32_t subset02[] = { 673,465,569,603,629,592,584,300,601,599,600,0 } ;
static int32_t subset03[] = { 90,85,83,84,65,87,76,46,0 } ;
static int32_t subset04[] = { 165,168,169,190,162,85,176,167,127,0 } ;
static int32_t subset05[] = { 224,275,278,249,277,279,289,295,139,0 } ;
static int32_t subset06[] = { 354,370,362,384,359,324,360,180,350,270,0 } ;
static int32_t subset07[] = { 599,595,557,298,448,596,577,667,597,588,602,0 } ;
static int32_t subset08[] = { 175,199,137,88,187,173,168,171,174,0 } ;
static int32_t subset09[] = { 93,187,196,144,185,178,186,202,182,0 } ;
static int32_t subset10[] = { 157,155,81,158,119,176,152,167,159,0 } ;
static int32_t subset11[] = { 184,165,159,166,163,167,174,124,83,0 } ;
static int32_t subset12[] = { 1211,1212,1287,605,1208,1189,1060,1216,1243,1200,908,1210,0 } ;
static int32_t subset13[] = { 339,299,153,305,282,304,313,306,302,228,0 } ;
static int32_t subset14[] = { 94,104,63,112,80,84,93,96,0 } ;
static int32_t subset15[] = { 41,88,82,85,61,74,83,81,0 } ;
static int32_t subset16[] = { 90,67,84,83,82,97,86,41,0 } ;
static int32_t subset17[] = { 299,303,151,301,291,302,307,377,333,280,0 } ;
static int32_t subset18[] = { 55,40,48,44,25,42,41,0 } ;
static int32_t subset19[] = { 1038,1188,1255,1184,594,890,1173,1151,1186,1203,1187,1195,0 } ;
static int32_t subset20[] = { 76,132,133,144,135,99,128,154,0 } ;
static int32_t subset21[] = { 77,46,108,81,85,84,93,83,0 } ;
static int32_t subset22[] = { 624,596,391,605,529,610,607,568,604,603,453,0 } ;
static int32_t subset23[] = { 83,167,166,189,163,174,160,165,133,0 } ;
static int32_t subset24[] = { 308,281,389,292,346,303,302,304,300,173,0 } ;
static int32_t subset25[] = { 593,1151,1187,1184,890,1040,1173,1186,1195,1255,1188,1203,0 } ;
static int32_t subset26[] = { 68,46,64,33,60,58,65,0 } ;
static int32_t subset27[] = { 65,43,88,87,86,99,93,90,0 } ;
static int32_t subset28[] = { 83,78,107,48,84,87,96,85,0 } ;
static int32_t subset29[] = { 1188,1173,1256,1038,1187,1151,890,1186,1184,1203,594,1195,0 } ;
static int32_t subset30[] = { 302,324,280,296,294,160,367,298,264,299,0 } ;
static int32_t subset31[] = { 521,760,682,687,646,664,342,698,692,686,672,0 } ;
static int32_t subset32[] = { 56,95,86,97,96,89,108,120,0 } ;
static int32_t subset33[] = { 344,356,262,343,340,382,337,175,361,330,0 } ;
static int32_t subset34[] = { 47,44,42,27,41,40,37,0 } ;
static int32_t subset35[] = { 139,155,161,158,118,166,154,156,78,0 } ;
static int32_t subset36[] = { 118,157,164,158,161,79,139,150,159,0 } ;
static int32_t subset37[] = { 299,292,371,150,300,301,281,303,306,262,0 } ;
static int32_t subset38[] = { 85,77,86,84,44,88,91,67,0 } ;
static int32_t subset39[] = { 88,85,84,44,65,91,76,86,0 } ;
static int32_t subset40[] = { 138,141,127,96,136,154,135,76,0 } ;
static int32_t subset41[] = { 292,308,302,346,300,324,304,305,238,166,0 } ;
static int32_t subset42[] = { 354,342,341,257,348,343,345,321,170,301,0 } ;
static int32_t subset43[] = { 84,178,168,167,131,170,193,166,162,0 } ;
static int32_t subset44[] = { 686,701,706,673,694,687,652,343,683,606,518,0 } ;
static int32_t subset45[] = { 295,293,301,367,296,279,297,263,323,159,0 } ;
static int32_t subset46[] = { 1038,1184,593,890,1188,1173,1187,1186,1195,1150,1203,1255,0 } ;
static int32_t subset47[] = { 343,364,388,402,191,383,382,385,288,374,0 } ;
static int32_t subset48[] = { 1187,1036,1183,591,1184,1175,888,1197,1182,1219,1115,1167,0 } ;
static int32_t subset49[] = { 151,291,307,303,345,238,299,323,301,302,0 } ;
static int32_t subset50[] = { 140,151,143,138,99,69,131,137,0 } ;
static int32_t subset51[] = { 29,44,42,59,41,36,40,0 } ;
static int32_t subset52[] = { 348,329,343,344,338,315,169,359,375,271,0 } ;
static int32_t subset53[] = { 48,39,34,37,50,40,41,0 } ;
static int32_t subset54[] = { 593,445,595,558,662,602,591,297,610,580,594,0 } ;
static int32_t subset55[] = { 686,651,681,342,541,687,691,707,604,675,699,0 } ;
static int32_t subset56[] = { 180,99,189,166,194,188,144,187,199,0 } ;
static int32_t subset57[] = { 321,349,335,343,377,176,265,356,344,332,0 } ;
static int32_t subset58[] = { 1151,1255,1195,1173,1184,1186,1188,1187,1203,593,1038,891,0 } ;
static int32_t subset59[] = { 90,88,100,83,62,113,80,89,0 } ;
static int32_t subset60[] = { 308,303,238,300,151,304,324,293,346,302,0 } ;
static int32_t subset61[] = { 59,38,50,41,42,35,40,0 } ;
static int32_t subset62[] = { 352,366,174,355,344,265,343,310,338,331,0 } ;
static int32_t subset63[] = { 91,89,93,90,117,85,60,106,0 } ;
static int32_t subset64[] = { 146,186,166,175,202,92,184,183,189,0 } ;
static int32_t subset65[] = { 82,67,96,44,80,79,88,76,0 } ;
static int32_t subset66[] = { 54,50,58,66,31,61,64,0 } ;
static int32_t subset67[] = { 343,266,344,172,308,336,364,350,359,333,0 } ;
static int32_t subset68[] = { 88,49,87,82,90,98,86,115,0 } ;
static int32_t subset69[] = { 20,47,49,51,54,48,40,0 } ;
static int32_t subset70[] = { 159,79,177,158,157,152,155,167,118,0 } ;
static int32_t subset71[] = { 1219,1183,1182,1115,1035,1186,591,1197,1167,887,1184,1175,0 } ;
static int32_t subset72[] = { 611,518,693,343,704,667,686,682,677,687,725,0 } ;
static int32_t subset73[] = { 607,599,634,305,677,604,603,580,452,605,591,0 } ;
static int32_t subset74[] = { 682,686,635,675,692,730,687,342,517,658,695,0 } ;
static int32_t subset75[] = { 662,296,573,598,592,584,553,593,595,443,591,0 } ;
static int32_t subset76[] = { 180,185,186,199,187,210,93,177,149,0 } ;
static int32_t subset77[] = { 197,136,179,185,156,182,180,178,99,0 } ;
static int32_t subset78[] = { 271,298,218,279,285,282,280,238,140,0 } ;
static int32_t subset79[] = { 1187,1151,890,593,1194,1188,1184,1173,1038,1186,1255,1203,0 } ;
static int32_t subset80[] = { 169,161,177,192,130,165,84,167,168,0 } ;
static int32_t subset81[] = { 50,42,43,41,66,39,36,0 } ;
static int32_t subset82[] = { 590,669,604,579,448,599,560,299,601,597,598,0 } ;
static int32_t subset83[] = { 174,191,206,179,184,142,177,180,90,0 } ;
static int32_t subset84[] = { 298,299,297,306,164,285,374,269,329,295,0 } ;
static int32_t subset85[] = { 181,172,162,138,170,195,86,169,168,0 } ;
static int32_t subset86[] = { 1184,1197,591,1182,1186,889,1167,1219,1183,1033,1115,1175,0 } ;
static int32_t subset87[] = { 644,695,691,679,667,687,340,681,770,686,517,0 } ;
static int32_t subset88[] = { 606,524,592,576,628,593,591,584,296,444,595,0 } ;
static int32_t subset89[] = { 94,127,154,138,135,74,136,141,0 } ;
static int32_t subset90[] = { 179,168,172,178,177,89,198,186,137,0 } ;
static int32_t subset91[] = { 302,299,291,300,298,149,260,305,280,370,0 } ;
static int32_t subset92[] = { 678,517,670,686,682,768,687,648,342,692,702,0 } ;
static int32_t subset93[] = { 302,290,304,376,333,303,306,298,279,153,0 } ;
static int32_t subset94[] = { 95,102,109,54,96,75,85,97,0 } ;
static int32_t subset95[] = { 150,154,146,78,152,151,162,173,119,0 } ;
static int32_t subset96[] = { 150,143,157,152,184,112,154,151,132,0 } ;
static int32_t subset97[] = { 36,41,54,40,25,44,42,0 } ;
static int32_t subset98[] = { 37,48,34,59,39,41,40,0 } ;
static int32_t subset99[] = { 681,603,638,611,584,303,454,607,606,605,596,0 } ;


static int32_t * tab105[101] = {
     subset00    ,subset01    ,subset02    ,subset03    ,subset04    ,subset05    ,subset06    ,subset07    ,subset08   ,subset09
    ,subset10    ,subset11    ,subset12    ,subset13    ,subset14    ,subset15    ,subset16    ,subset17    ,subset18    ,subset19
    ,subset20    ,subset21    ,subset22    ,subset23    ,subset24    ,subset25    ,subset26    ,subset27    ,subset28    ,subset29
    ,subset30    ,subset31    ,subset32    ,subset33    ,subset34    ,subset35    ,subset36    ,subset37    ,subset38    ,subset39
    ,subset40    ,subset41    ,subset42    ,subset43    ,subset44    ,subset45    ,subset46    ,subset47    ,subset48    ,subset49
    ,subset50    ,subset51    ,subset52    ,subset53    ,subset54    ,subset55    ,subset56    ,subset57    ,subset58    ,subset59
    ,subset60    ,subset61    ,subset62    ,subset63    ,subset64    ,subset65    ,subset66    ,subset67    ,subset68    ,subset69
    ,subset70    ,subset71    ,subset72    ,subset73    ,subset74    ,subset75    ,subset76    ,subset77    ,subset78    ,subset79
    ,subset80    ,subset81    ,subset82    ,subset83    ,subset84    ,subset85    ,subset86    ,subset87    ,subset88    ,subset89
    ,subset90    ,subset91    ,subset92    ,subset93    ,subset94    ,subset95    ,subset96    ,subset97    ,subset98    ,subset99
    ,NULL
} ;

int32_t * * P105_GetData(void) {
    return tab105 ;
}
