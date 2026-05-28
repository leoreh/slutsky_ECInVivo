// Functions:
// (1) initPins - initialize pins. Called once during setup.
// (2) sensThr - determine sensors threshold. Called once during setup.
// (3) readSns - Read sensors value and return if threshold was crossed. Called constantely from detectRun.
// (4) giveReward - opens relavent solenoid. Called from detectRun.
// (5) printSens - prints sensor readings on screen. Used only for debugging.
// (6) detectRun - Calls readSns. If thresholds are crossed in a specific sequence than calls giveReward. Called from Arduino's loop.

// Fixed parameters:
const int IRLED_pins[] = {2};
const int sens_pins[] = {A0, A1, A2, A3};
const int sol_pins[] = {8, 9};
const int thr_pins[] = {3, 4, 5, 6};   // pins that report to DAQ if a sensor has been crossed
const int man_sw[] = {7};
const int numSens = sizeof(sens_pins) / sizeof(int);;
const int numSol = sizeof(sol_pins) / sizeof(int);;
const int numThrPins = sizeof(thr_pins) / sizeof(int);;
int thr[] = {0, 0, 0, 0};
int r2l[] = {0, 0, 0};
int l2r[] = {0, 0, 0};
int openDuration  = 30;
int runNum = 0;

void setup() {
  Serial.begin(9600);
  initPins();
  sensorThr();
}

void loop() {
  detectRun();
}

void initPins() {
  pinMode(IRLED_pins[0], OUTPUT);
  pinMode(man_sw[0], INPUT_PULLUP);
  digitalWrite(IRLED_pins[0], HIGH);
  for (int i = 0; i < numSens; i++) {
    pinMode(sens_pins[i], INPUT);
  }
  for (int i = 0; i < numSol; i++) {
    pinMode(sol_pins[i], OUTPUT);
  }
  for (int i = 0; i < numThrPins; i++) {
    pinMode(thr_pins[i], OUTPUT);
  }
}

void sensorThr() {
  delay(500);
  for (int i = 0; i < numSens; i++) {
    thr[i] = analogRead(sens_pins[i]) + 35;
    Serial.print(i + 1);
    Serial.print(" :");
    Serial.println(thr[i]);
  }
}

int readSens() {
  for (int i = 0; i < numSens; i++) {
    int value = analogRead(sens_pins[i]);
    delay(2);
    // Serial.print(i+1);
    // Serial.print(" :");
    // Serial.println(value);
    if (value > thr[i]) {
      digitalWrite(thr_pins[i], HIGH);
      return i;
    }
    else {
      digitalWrite(thr_pins[i], LOW);
      return -1;
    }
  }
}

void giveReward(int solNo) {
  digitalWrite(sol_pins[solNo], HIGH);
  delay(openDuration);
  digitalWrite(sol_pins[solNo], LOW);
  runNum ++;
  r2l[0] = 0;
  r2l[1] = 0;
  r2l[2] = 0;
  l2r[0] = 0;
  l2r[1] = 0;
  l2r[2] = 0;
}

void printSens() {
  for (int i = 0; i < numSens; i++) {
    Serial.print(i + 1);
    Serial.print(" :");
    Serial.println(analogRead(sens_pins[i]));
    delay(1000);
  }
}

void detectRun() {
  // right to left
  if (readSens() == 0)
    r2l[0] = 1;
  if (readSens() == 1 && r2l[0] == 1)
    r2l[1] = 1;
  if (readSens() == 2 && r2l[1] == 1)
    r2l[2] = 1;
  if (readSens() == 3 && r2l[2] == 1)
    giveReward(1);
  // left to right
  if (readSens() == 3)
    l2r[0] = 1;
  if (readSens() == 2 && l2r[0] == 1)
    l2r[1] = 1;
  if (readSens() == 1 && l2r[1] == 1)
    l2r[2] = 1;
  if (readSens() == 0 && l2r[2] == 1)
    giveReward(0);
}





