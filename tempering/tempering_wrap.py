def init(infile, outfile):
    return _tempering.tempering_init_from_plugin(infile, len(infile), outfile, len(outfile))
def update(iteration, energy, temperature):
    ierror, new_temperature = _tempering.tempering_dyna_from_plugin(iteration, energy, temperature)
    return new_temperature;
def done():
    return _tempering.tempering_done_from_plugin()
