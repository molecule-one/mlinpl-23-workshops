# Workshop scope

You are drug hunter. During the workshops, we will code various version of an active learning loop that orders experiments in the lab.

Only easy to synthesize are accepted! You have limited # of calls but can synthesize any number

## Agenda

1. 9:00 - 9:20 Introduction 
2. 9:30 - 10:30 Work on task 1,2,3 
5. 10:30 - 11:00: What did we learn? Discussion
6. 11:00 - 11:15 Break

## 1. Tasks

Your goal will be to create various implementation of active learning. Technically, you will be tasked with code various implementations of the `Loop` base class.

### 1.1. Seed step

Submit first step of active learning to the server and see it on the leaderboard
Given that we have no prior knowledge about the target, we will start with random search.

### 1.2. Exploration step

Implement mutation based active learning.

Tune hyperparameters to achieve 95% activity on the DRD2 protein in 1000 samples.

Answer: What makes it work, what doesn't? Short discussion will follow.

### 1.3. Compare random and mutation

Implement comparison script/function

### 1.4. ML guided active learning

(TODO)

Investigate synthesizability random and mutation for 500 budget compounds using mol2grid. What do you observe?
What is the largest

## 2. Working with the server

Assumes you have a running instance accessible at `127.0.0.1:5000`.

### 2.1. Starting the server

```commandline
PYTHONPATH=$PYTHONPATH:`pwd` python server/start.py
```

### 2.2. Reseting

CURL:
```commandline
curl -X POST -H "Content-Type: application/json" -d '{}' http://127.0.0.1:5000/reset
```

## 3. Running solution

### 3.1 Run random search on server

```commandline
 python solutions/run.py -b 1000 -w random -t DRD2_server -u test-0
```
