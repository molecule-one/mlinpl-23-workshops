# Workshop scope

You are drug hunter. During the workshops, we will code various version of an active learning loop that orders experiments in the lab.

Only easy to synthesize are accepted! You have limited # of calls but can synthesize any number

## 1. Tasks

Your goal will be to create various implementation of active learning. Technically, you will be tasked with code various implementations of the `Loop` base class.

### 1.1. Seed step

Submit first step of active learning to the server and see it on the leaderboard
Given that we have no prior knowledge about the target, we will start with random search.

### 1.2. Exploration step

Implement mutation based active learning.

### 1.3. ML guided active learning

(TODO)

Investigate synthesizability random and mutation for 500 budget compounds using mol2grid. What do you observe?
What is the largest

## 2. Working with the server

Assumes you have a running instance accessible at `localhost:5000`.

### 2.1. Starting the server

### 2.2. Reseting

CURL:
```commandline
curl -X POST -H "Content-Type: application/json" -d '{}' http://127.0.0.1:5000/reset
```

## 2.3. Adding results

Python:

CURL:

```commandline
curl -X POST -H "Content-Type: application/json" \
-d '{
    "token": "unique_token_1",
    "metrics": {
        "metric1": 10,
        "metric2": 20
    }
}' \
http://127.0.0.1:5000/add_result

```