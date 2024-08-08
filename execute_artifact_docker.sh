#!/bin/bash

REPO_ROOT=$(git rev-parse --show-toplevel)
cd ${REPO_ROOT}

echo '## Acquiring docker image'
docker pull ucdavisplse/prose-env:sc24

if ! ls paper_results > /dev/null 2>&1
then
    cmd="source .venv/bin/activate && \
        echo '## Moving paper data from the docker container to ${REPO_ROOT}/paper_results' && \
        cd /root/artifact && \
        pv -f /root/results-*.tar.gz | tar -xz && \
        mv results-* paper_results && \        echo '## Generating interactive plots from paper data and saving to ${REPO_ROOT}' && \
        cd /root/artifact/paper_results && \
        for EXPERIMENT_NAME in *; do \
            SCRIPT_PATH=\$(cd /root/artifact && find ./submodules/\${EXPERIMENT_NAME} -name 'generate_figures.py') && \
            [ -n \${SCRIPT_PATH} ] && \
                ln -sf /root/artifact/\${SCRIPT_PATH} /root/artifact/paper_results/\${EXPERIMENT_NAME}/generate_figures.py && \
                pushd \${EXPERIMENT_NAME} > /dev/null 2>&1 && \
                echo '         executing '\${SCRIPT_PATH}' in ./paper_results/'\${EXPERIMENT_NAME} && \
                python3 generate_figures.py > /dev/null 2>&1 && \
                rm generate_figures.py && \
                mv *.html /root/artifact && \
                popd > /dev/null 2>&1 ; \
        done"

    docker run --rm -t -v ${REPO_ROOT}:/root/artifact ucdavisplse/prose-env:sc24 bash -c "${cmd}"
fi

cmd="source .venv/bin/activate && \
    echo && \
    echo '## Applying tool to funarc example' && \
    echo && \
    cd /root/artifact/submodules/funarc && \
    make reset > /dev/null 2>&1 && \
    prose_search.py -s setup.ini && \
    python3 generate_figures.py && \
    mv *.html /root/artifact"

docker run --rm -t -v ${REPO_ROOT}:/root/artifact ucdavisplse/prose-env:sc24 bash -c "${cmd}"