[source dir]      [path_inputs]
    +                  +
    |                  |
    |                  |
    v                  v               +Experiment1 of neuron1
 Archon  <--------+ Chorus <-----------+Experiment1 of neuron2
                      +                +Experiment1 of neuron3
                      |                           ^        ...
                      |                           |
                      |                           |                          +Piece1       +Sound
                      |                           +--------------------------+Piece2 <-----+Spike
                      |                           |                          +Piece3       +Trigger
                      |                           |                              +
                      |                           |                              |
                      +                           v                              |
                                       +Experiment1                              |
      Sultan <---------+Neuron <-------+Experiment2 of one neuron                |
                          +            +Experiment3                              |
                      +   |             +       ...                              |
                      |   |             |                                        |
                      |   |             |                                        |
                      |   |             |                                        |
                      |   |             |                                        |
                      |   |             v                                        |
                      |   +----> EphysAnalysis  <--------------------------------+
                      |                ^
                      |                |         <inherit>
                      |                |
                      |                |
                      +----------------+
