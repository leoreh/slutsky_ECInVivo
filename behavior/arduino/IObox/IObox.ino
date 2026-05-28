/*
  IO box controlling open ephys digital signals
  LH 07 jun 20

  Functions:
  (1) initPins      initialize pins. Called once during setup.
  (2) cam           send pulse to camera (frame acquisition). called in loop.
  (3) stim          send pulse to stimulator according to pre-defined protocols:
                    1 - IO; single pulse repeated every isi sec {10 = 0.1 Hz}
                    2 - STP; 5 pulses of 50 Hz repeated once every 15 seconds
                    3 - LTP; five bursts of 50 x 200 Hz stimulations once a second, repeated six times once a minute
*/

// Fixed parameters:
const float camF = 10;              // camera frequency [fps]
const int pulseDur = 1;             // output pulse duration [ms]
const int cam_pins[] = {8, 11};
const int stim_pins[] = {9, 12};
const int camSwch_pin = 2;
const int stimSwch_pin = 3;
const int abortSwch_pin = 7;        // force stops ltp stimulation
int rep[3], i, ii, iii;

// custom params
int stimProtocol = 2;               // see above
int overrideDef = 0;                // use stimF and isi defined above instead of default values
float stimF = 50;                   // to override stp
float isi = 2;                      // time between stims [s]. for io and stp only

void setup() {
  Serial.begin(9600);
  initPins();
}

void loop() {
  int camSwch_val = digitalRead(camSwch_pin);
  int stimSwch_val = digitalRead(stimSwch_pin);
  if (camSwch_val < 1) {
    cam();
  }
  if (stimSwch_val < 1) {
    stim();
  }
}

void initPins() {
  pinMode(camSwch_pin, INPUT_PULLUP);
  pinMode(stimSwch_pin, INPUT_PULLUP);
  pinMode(abortSwch_pin, INPUT_PULLUP);
  for (i = 0; i < 2; i++) {
    pinMode(cam_pins[i], OUTPUT);
  }
  for (i = 0; i < 2; i++) {
    pinMode(stim_pins[i], OUTPUT);
  }
}

void cam() {
  digitalWrite(cam_pins[0], HIGH);
  digitalWrite(cam_pins[1], HIGH);
  delay(pulseDur);
  digitalWrite(cam_pins[0], LOW);
  digitalWrite(cam_pins[1], LOW);
  delay(1 / camF * 1000 - pulseDur);
}

void stim() {
  switch (stimProtocol) {
    case 1:     // io
      if (overrideDef < 1) {
        isi = 20;;
      }
      stimF = 0.1;
      digitalWrite(stim_pins[0], HIGH);
      digitalWrite(stim_pins[1], HIGH);
      delay(pulseDur);
      digitalWrite(stim_pins[0], LOW);
      digitalWrite(stim_pins[1], LOW);
      delay(isi * 1000 - pulseDur);
      break;

    case 2:      // stp
      rep[0] = 3;
      if (overrideDef < 1) {
        stimF = 50;
        isi = 30;;
      }
      for (i = 0; i < rep[0]; i++) {
        digitalWrite(stim_pins[0], HIGH);
        digitalWrite(stim_pins[1], HIGH);
        delay(pulseDur);
        digitalWrite(stim_pins[0], LOW);
        digitalWrite(stim_pins[1], LOW);
        delay(1 / stimF * 1000 - pulseDur);
      }
      delay(isi  * 1000);
      break;

    case 3:      // ltp
      int rep[] = {50, 5, 6};
      stimF = 200;
      for (i = 0; i < rep[2]; i++) {
        for (ii = 0; ii < rep[1]; ii++) {
          for (iii = 0; iii < rep[0]; iii++) {
            digitalWrite(stim_pins[0], HIGH);
            digitalWrite(stim_pins[1], HIGH);
            delay(pulseDur);
            digitalWrite(stim_pins[0], LOW);
            digitalWrite(stim_pins[1], LOW);
            delay(1 / stimF * 1000 - pulseDur);
          }
          delay(1000);
          int abortSwch_val = digitalRead(abortSwch_pin);
          if (abortSwch_val < 1) {
            return;
          }
        }
        delay(60000);
      }
      break;

    default:
      break;
  }
}
