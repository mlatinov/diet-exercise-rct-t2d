library(dagitty)

dag <- dagitty('
dag {
bb="-10,-10,10,10"

PIP [exposure, pos="-7,0"]

Age [pos="-9,-7"]
Sex [pos="-7,-7"]
Marital [pos="-5,-7"]
Height [pos="-3,-7"]

Qpre [pos="-5,-3"]
SDSCAPre [pos="-3,-3"]
DietPre [pos="-1,-3"]
CarbPre [pos="1,-3"]
EnergyPre [pos="3,-3"]
ExPre [pos="5,-3"]
ExIntPre [pos="6,-4"]
ExDurPre [pos="7,-3"]

WtPre [pos="-3,-5"]
BMIPre [pos="-1,-5"]
WaistPre [pos="1,-5"]
HipPre [pos="3,-5"]

BarTime [pos="-9,2"]
BarCount [pos="-7,2"]
BarMot [pos="-5,2"]
BarSocial [pos="-3,2"]
BarFood [pos="-1,2"]
BarHome [pos="1,2"]

Qpost [pos="-5,4"]
SDSCAPost [pos="-3,4"]
Adhere [pos="-1,4"]
ExAdh [pos="1,4"]
DietPost [pos="3,4"]
CarbPost [pos="5,4"]
EnergyPost [pos="7,4"]
ExPost [pos="9,4"]
ExIntPost [pos="9,5"]
ExDurPost [pos="9,6"]

WtPost [outcome, pos="-3,7"]
BMIPost [outcome, pos="-1,7"]
WaistPost [outcome, pos="1,7"]
HipPost [outcome, pos="3,7"]

Age -> CarbPre
Age -> DietPre
Age -> EnergyPre
Age -> ExDurPre
Age -> ExIntPre
Age -> ExPre
Age -> HipPre
Age -> Qpre
Age -> SDSCAPre
Age -> WaistPre
Age -> WtPre
Age -> BMIPre

Sex -> CarbPre
Sex -> EnergyPre
Sex -> ExIntPre
Sex -> ExPre
Sex -> HipPre
Sex -> WaistPre
Sex -> WtPre

Marital -> CarbPre
Marital -> DietPre
Marital -> EnergyPre
Marital -> ExPre
Marital -> SDSCAPre

Height -> BMIPost
Height -> BMIPre
Height -> WtPost
Height -> WtPre

WtPre -> BMIPre
WtPost -> BMIPost

CarbPre -> EnergyPre
SDSCAPre -> CarbPre
SDSCAPre -> DietPre
SDSCAPre -> ExPre
SDSCAPre -> WaistPre
SDSCAPre -> WtPre
ExIntPre -> ExPre
ExDurPre -> ExPre
Qpre -> DietPre
Qpre -> SDSCAPre
CarbPre -> WtPre
EnergyPre -> WtPre
DietPre -> WtPre
ExPre -> WtPre

WtPre -> WtPost
BMIPre -> BMIPost
WaistPre -> WaistPost
HipPre -> HipPost
CarbPre -> CarbPost
EnergyPre -> EnergyPost
DietPre -> DietPost
ExPre -> ExPost
ExIntPre -> ExIntPost
ExDurPre -> ExDurPost
SDSCAPre -> SDSCAPost
Qpre -> Qpost

Age -> BarTime
Age -> BarHome
Sex -> BarTime
Sex -> BarHome
Marital -> BarHome
Marital -> BarSocial

PIP -> BarTime
PIP -> BarCount
PIP -> BarMot
PIP -> BarSocial
PIP -> BarFood
PIP -> BarHome

Qpre -> BarCount
SDSCAPre -> BarMot

BarTime -> CarbPost
BarTime -> DietPost
BarTime -> ExPost
BarTime -> Adhere
BarCount -> CarbPost
BarCount -> DietPost
BarCount -> Adhere
BarMot -> SDSCAPost
BarMot -> Adhere
BarMot -> ExAdh
BarMot -> ExPost
BarSocial -> CarbPost
BarSocial -> EnergyPost
BarSocial -> Adhere
BarFood -> CarbPost
BarFood -> EnergyPost
BarFood -> DietPost
BarFood -> Adhere
BarHome -> ExPost
BarHome -> ExAdh
BarHome -> Adhere
BarHome -> SDSCAPost

PIP -> Adhere
PIP -> CarbPost
PIP -> DietPost
PIP -> EnergyPost
PIP -> ExAdh
PIP -> ExDurPost
PIP -> ExIntPost
PIP -> ExPost
PIP -> Qpost
PIP -> SDSCAPost

Qpost -> Adhere
Qpost -> DietPost
Qpost -> SDSCAPost
SDSCAPost -> Adhere
SDSCAPost -> CarbPost
SDSCAPost -> DietPost
SDSCAPost -> ExAdh
SDSCAPost -> ExPost
Adhere -> CarbPost
Adhere -> DietPost
Adhere -> EnergyPost
ExAdh -> ExDurPost
ExAdh -> ExIntPost
ExAdh -> ExPost
CarbPost -> EnergyPost
ExIntPost -> ExPost
ExDurPost -> ExPost

Adhere -> WaistPost
Adhere -> WtPost
CarbPost -> WtPost
DietPost -> WtPost
EnergyPost -> WtPost
ExPost -> WaistPost
ExPost -> WtPost
SDSCAPost -> WaistPost
SDSCAPost -> WtPost

WtPost -> HipPost
WtPost -> WaistPost
WaistPost -> HipPost
}
')
adjustmentSets(dag, exposure = "PIP", outcome = "HipPre", effect = "total")
adjustmentSets(dag, exposure = "PIP", outcome = "HipPre", effect = "direct")
