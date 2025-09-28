# Rate Limit Guide

## Why am I seeing rate limit errors?

The Fastalysis app uses AI models from Groq and other providers that have usage limits to prevent abuse. When you make too many requests in a short period, you may encounter rate limiting.

## How to fix rate limit issues:

1. **Wait a few minutes before trying again**
   * Rate limits typically reset after a short cooling-off period (2-5 minutes)

2. **Try a different model**
   * Smaller models like `llama-3.1-8b-instant` have higher rate limits
   * BioGPT uses a different API provider and may be available when Groq is rate-limited

3. **Keep your queries concise and specific**
   * Long conversations with many back-and-forth messages are more likely to hit rate limits

4. **Use the built-in analysis tools**
   * The BLAST search and sequence analysis tools don't use the rate-limited API
   * Try using these specialized tools instead of asking the chat for analysis

5. **If problems persist:**
   * Try again during off-peak hours
   * Check if there are any API service outages
   * Make sure your API keys are configured correctly in `.env.local`

## Technical details

Rate limits are typically based on:
* Tokens per minute
* Requests per minute
* Daily usage quotas

Different models have different limits, so switching models can help when one is rate-limited.