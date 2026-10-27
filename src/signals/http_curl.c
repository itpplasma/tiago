#define _POSIX_C_SOURCE 200809L
#include <curl/curl.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* C owns varargs and platform clock/locking ABI; Fortran owns the client. */
static pthread_once_t curl_once = PTHREAD_ONCE_INIT;
static int curl_status;
static void initialize_curl(void) { curl_status = curl_global_init(CURL_GLOBAL_DEFAULT); }
int tiago_curl_global(void) {
    pthread_once(&curl_once, initialize_curl);
    return curl_status;
}
double tiago_http_clock(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}
void tiago_http_sleep(double seconds) {
    if (seconds <= 0) return;
    struct timespec ts = {(time_t)seconds, (long)((seconds - (time_t)seconds)*1e9)};
    while (nanosleep(&ts, &ts) && errno == EINTR) {}
}
typedef struct host_gate {
    char *host;
    pthread_mutex_t mutex;
    double last_start, interval, embargo;
    struct host_gate *next;
} host_gate;
static pthread_mutex_t hosts_mutex = PTHREAD_MUTEX_INITIALIZER;
static host_gate *hosts;

void *tiago_http_acquire(const char *url, double interval, double deadline) {
    CURLU *u = curl_url();
    char *host = NULL, *scheme = NULL, *user = NULL;
    if (!u) return NULL;
    int valid = !curl_url_set(u, CURLUPART_URL, url, 0) &&
        !curl_url_get(u, CURLUPART_SCHEME, &scheme, 0) &&
        !curl_url_get(u, CURLUPART_HOST, &host, 0);
    if (valid) valid = !strcmp(scheme, "http") || !strcmp(scheme, "https");
    /* Credentials belong to an explicit private transport, never a URL. */
    if (!curl_url_get(u, CURLUPART_USER, &user, 0)) valid = 0;
    curl_free(user);
    curl_free(scheme);
    curl_url_cleanup(u);
    if (!valid) { curl_free(host); return NULL; }
    for (char *c = host; *c; ++c) *c = (char)tolower((unsigned char)*c);
    pthread_mutex_lock(&hosts_mutex);
    host_gate *g = hosts;
    while (g && strcmp(g->host, host)) g = g->next;
    if (!g) {
        g = calloc(1, sizeof(*g));
        if (g) {
            g->host = strdup(host);
            if (!g->host) { free(g); g = NULL; }
        }
        if (g) {
            pthread_mutex_init(&g->mutex, NULL);
            g->next = hosts;
            hosts = g;
        }
    }
    pthread_mutex_unlock(&hosts_mutex);
    curl_free(host);
    if (!g) return NULL;
    while (pthread_mutex_trylock(&g->mutex)) {
        if (tiago_http_clock() >= deadline) return NULL;
        tiago_http_sleep(0.001);
    }
    double spacing = interval > g->interval ? interval : g->interval;
    double ready = g->last_start + spacing;
    if (g->embargo > ready) ready = g->embargo;
    double now = tiago_http_clock();
    if (ready >= deadline || now >= deadline) {
        pthread_mutex_unlock(&g->mutex);
        return NULL;
    }
    tiago_http_sleep(ready - now);
    g->last_start = tiago_http_clock();
    g->interval = interval;
    return g;
}
void tiago_http_release(void *gate, double delay) {
    host_gate *g = gate;
    if (delay > 0) g->embargo = tiago_http_clock() + delay;
    pthread_mutex_unlock(&g->mutex);
}
int tiago_curl_setup(void *handle, const char *url, const char *agent,
                     long timeout, long connect_timeout,
                     size_t (*writer)(char *, size_t, size_t, void *), void *data) {
    CURL *c = handle;
    CURLcode rc;
#define SET(opt, val) do { rc = curl_easy_setopt(c, opt, val); if (rc) return rc; } while (0)
    SET(CURLOPT_URL, url);
    SET(CURLOPT_USERAGENT, agent);
    SET(CURLOPT_TIMEOUT_MS, timeout);
    SET(CURLOPT_CONNECTTIMEOUT_MS, connect_timeout);
    SET(CURLOPT_NOSIGNAL, 1L);
    SET(CURLOPT_SSL_VERIFYPEER, 1L);
    SET(CURLOPT_SSL_VERIFYHOST, 2L);
    SET(CURLOPT_FOLLOWLOCATION, 0L);
#if LIBCURL_VERSION_NUM >= 0x075500
    SET(CURLOPT_PROTOCOLS_STR, "http,https");
#else
    SET(CURLOPT_PROTOCOLS, (long)(CURLPROTO_HTTP | CURLPROTO_HTTPS));
#endif
    SET(CURLOPT_WRITEFUNCTION, writer);
    SET(CURLOPT_WRITEDATA, data);
#undef SET
    return CURLE_OK;
}
int tiago_curl_info(void *handle, long *code, double *retry) {
    curl_off_t seconds = 0;
    CURLcode rc = curl_easy_getinfo(handle, CURLINFO_RESPONSE_CODE, code);
    if (rc) return rc;
    rc = curl_easy_getinfo(handle, CURLINFO_RETRY_AFTER, &seconds);
    *retry = seconds > 0 ? (double)seconds : 0;
    return rc;
}
