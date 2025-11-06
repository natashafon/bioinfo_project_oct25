# Start from Ubuntu
FROM ubuntu:22.04

# Avoid prompts during install
ENV DEBIAN_FRONTEND=noninteractive

# Install R and dependencies
RUN apt-get update && apt-get install -y \
    r-base \
    gdebi-core \
    wget \
    sudo \
    && apt-get clean

# Install RStudio Server
RUN wget https://download2.rstudio.org/server/jammy/amd64/rstudio-server-2024.04.2-764-amd64.deb -O /tmp/rstudio-server.deb && \
    gdebi -n /tmp/rstudio-server.deb && \
    rm /tmp/rstudio-server.deb

# Create RStudio user and password
RUN useradd -m rstudio && echo "rstudio:rstudio" | chpasswd && adduser rstudio sudo

# Install essential R packages
RUN R -q -e "install.packages(c('devtools'), repos='https://cloud.r-project.org')"
RUN R -q -e "install.packages(c('data.table', 'dplyr', 'ggplot2'), repos='https://cloud.r-project.org')"

# Set working directory
WORKDIR /home/project_oct

# Expose RStudio Server port
EXPOSE 8787

# Default command: start RStudio Server
CMD ["/usr/lib/rstudio-server/bin/rserver", "--server-daemonize=0"]
