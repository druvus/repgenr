// The run-level meta map: one per pipeline run, carried on every channel.
//
// id   -- a slug of the selection target, taken from the argument strings the
//         user passes to the metadata stages, so published work directories
//         and task tags name the taxon. Falls back to the mode.
//         Bacterial: the value after `-ts` in metadata_args, else after
//         `-tg`, else after `-tf` (most specific rank first). Viral: the
//         value after `-tg` in vgenome_args, else after `--target` in
//         vmetadata_args. Reads: as bacterial, over reads_args.
//         `target_after(args, flag)` returns the token after `flag` in
//         `args`, or null when `flag` is absent.
// mode -- 'bacterial', 'viral' or 'reads'.

def target_after(String args, String flag) {
    def tokens = args ? args.tokenize() : []
    def i = tokens.indexOf(flag)
    return (i >= 0 && i + 1 < tokens.size()) ? tokens[i + 1] : null
}

def run_meta(Map params) {
    def target = null
    if (params.mode == 'viral') {
        target = target_after(params.vgenome_args ?: '', '-tg')
            ?: target_after(params.vmetadata_args ?: '', '--target')
    }
    else {
        def args = (params.mode == 'reads' ? params.reads_args : params.metadata_args) ?: ''
        target = target_after(args, '-ts')
            ?: target_after(args, '-tg')
            ?: target_after(args, '-tf')
    }
    def stripped = (target ?: params.mode).replaceAll('^[\'"]+|[\'"]+$', '')
    def id = stripped.toLowerCase().replaceAll('[^a-z0-9]+', '_').replaceAll('^_+|_+$', '')
    return [id: id, mode: params.mode]
}
